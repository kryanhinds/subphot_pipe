"""
Quick seeing measurement from FITS files.

Measures PSF FWHM by detecting bright stars, fitting 2D Gaussians,
and converting to arcseconds.

Usage:
    from seeing_measure import measure_seeing_batch

    df['seeing'] = measure_seeing_batch(df['filename'], df['instrument'])
"""

import warnings
import numpy as np
from pathlib import Path
from typing import Union, List, Tuple, Optional, Dict
from astropy.io import fits
from scipy.ndimage import gaussian_filter
from scipy.optimize import curve_fit
import pandas as pd


# Pixel scales (arcsec/pixel) for different instruments
PIXEL_SCALES = {
    'NOT/ALFOSC': 0.189,
    'TJO/MEIA2': 0.304,
    'LT/IOO': 0.278,
    'SEDM': 1.01,
    'HCT': 0.3,
}


def get_pixel_scale(instrument: str) -> float:
    """
    Get pixel scale for an instrument.

    Parameters
    ----------
    instrument : str
        Instrument name (e.g., 'NOT/ALFOSC')

    Returns
    -------
    pixel_scale : float
        Arcsec per pixel
    """
    # Exact match
    if instrument in PIXEL_SCALES:
        return PIXEL_SCALES[instrument]

    # Partial match
    for instr, scale in PIXEL_SCALES.items():
        if instr.split('/')[0] in instrument or instr.split('/')[1] in instrument:
            return scale

    warnings.warn(f"Unknown instrument '{instrument}', using default 0.3 arcsec/pix")
    return 0.3


def gaussian_2d(coords: np.ndarray, amplitude: float, xo: float, yo: float,
                sigma_x: float, sigma_y: float, theta: float, offset: float) -> np.ndarray:
    """
    2D Gaussian function.

    Parameters
    ----------
    coords : ndarray
        Flattened (x, y) coordinate array
    amplitude : float
        Peak amplitude
    xo, yo : float
        Center coordinates
    sigma_x, sigma_y : float
        Standard deviations
    theta : float
        Rotation angle
    offset : float
        Background offset

    Returns
    -------
    result : ndarray
        Flattened Gaussian values
    """
    x, y = coords
    xo = float(xo)
    yo = float(yo)

    a = (np.cos(theta) ** 2) / (2 * sigma_x ** 2) + (np.sin(theta) ** 2) / (2 * sigma_y ** 2)
    b = -(np.sin(2 * theta)) / (4 * sigma_x ** 2) + (np.sin(2 * theta)) / (4 * sigma_y ** 2)
    c = (np.sin(theta) ** 2) / (2 * sigma_x ** 2) + (np.cos(theta) ** 2) / (2 * sigma_y ** 2)

    g = offset + amplitude * np.exp(-(a * (x - xo) ** 2 + 2 * b * (x - xo) * (y - yo) + c * (y - yo) ** 2))
    return g


def detect_stars(image: np.ndarray, threshold: float = 3.0,
                kernel_size: int = 5) -> List[Tuple[float, float, float]]:
    """
    Detect bright stars in image using simple thresholding.

    Parameters
    ----------
    image : ndarray
        2D image array
    threshold : float
        Detection threshold in sigma above background
    kernel_size : int
        Kernel size for background estimation

    Returns
    -------
    stars : list of (x, y, brightness)
        Detected star positions and peak values
    """
    # Remove NaNs and Infs
    image = np.nan_to_num(image, nan=0, posinf=0, neginf=0)

    # Estimate background
    blurred = gaussian_filter(image, sigma=kernel_size)
    background = np.median(blurred)
    noise = np.std(image - blurred)

    if noise == 0:
        warnings.warn("Image has zero noise, cannot detect stars")
        return []

    # Detect peaks above threshold
    detection_map = (image - background) / noise > threshold
    y_peaks, x_peaks = np.where(detection_map)

    if len(x_peaks) == 0:
        return []

    # Get peak values
    peak_values = image[y_peaks, x_peaks]

    # Cluster nearby peaks (simple approach: group within 10 pixels)
    stars = []
    used = set()

    for idx in np.argsort(-peak_values):  # Process brightest first
        if idx in used:
            continue

        x, y = x_peaks[idx], y_peaks[idx]
        brightness = peak_values[idx]

        # Mark nearby pixels as used
        mask = (np.abs(x_peaks - x) < 10) & (np.abs(y_peaks - y) < 10)
        used.update(np.where(mask)[0])

        stars.append((float(x), float(y), float(brightness)))

    return stars


def measure_psf_fwhm(image: np.ndarray, x: float, y: float,
                     box_size: int = 25) -> Optional[float]:
    """
    Measure PSF FWHM by fitting a 2D Gaussian to a star.

    Parameters
    ----------
    image : ndarray
        2D image array
    x, y : float
        Star center position
    box_size : int
        Size of cutout box around star (pixels)

    Returns
    -------
    fwhm : float or None
        FWHM in pixels, or None if fit failed
    """
    # Extract cutout
    x0, x1 = int(x - box_size // 2), int(x + box_size // 2)
    y0, y1 = int(y - box_size // 2), int(y + box_size // 2)

    # Bounds check
    if x0 < 0 or x1 >= image.shape[1] or y0 < 0 or y1 >= image.shape[0]:
        return None

    cutout = image[y0:y1, x0:x1]

    if cutout.size == 0:
        return None

    # Create coordinate grid
    yy, xx = np.meshgrid(np.arange(cutout.shape[0]), np.arange(cutout.shape[1]),
                         indexing='ij')
    coords = np.array([xx.ravel(), yy.ravel()])

    # Initial guess
    amplitude = np.max(cutout) - np.median(cutout)
    xo, yo = box_size // 2, box_size // 2
    sigma_init = 2.0
    offset = np.median(cutout)

    initial_guess = [amplitude, xo, yo, sigma_init, sigma_init, 0, offset]

    try:
        popt, _ = curve_fit(
            gaussian_2d, coords, cutout.ravel(),
            p0=initial_guess,
            maxfev=10000,
            ftol=1e-4,
        )

        # Extract sigma values and convert to FWHM
        sigma_x, sigma_y = popt[3], popt[4]
        sigma_mean = np.sqrt(sigma_x * sigma_y)
        fwhm = 2.355 * sigma_mean  # FWHM = 2.355 * sigma for Gaussian

        return fwhm

    except (RuntimeError, ValueError):
        return None


def measure_seeing_file(filename: str, instrument: str,
                       threshold: float = 3.0,
                       num_stars: int = 10,
                       verbose: bool = False) -> Optional[float]:
    """
    Measure seeing (FWHM) from a FITS file.

    Parameters
    ----------
    filename : str
        Path to FITS file
    instrument : str
        Instrument name for pixel scale lookup
    threshold : float
        Star detection threshold (sigma)
    num_stars : int
        Number of brightest stars to use for seeing estimate
    verbose : bool
        Print progress

    Returns
    -------
    seeing_arcsec : float or None
        Measured seeing in arcseconds, or None if measurement failed
    """
    try:
        with fits.open(filename) as hdul:
            image = hdul[0].data

        if image is None or image.size == 0:
            return None

        # Clean image
        image = np.nan_to_num(image, nan=0, posinf=0, neginf=0)

        # Detect stars
        stars = detect_stars(image, threshold=threshold)

        if len(stars) < 3:
            if verbose:
                print(f"  Warning: Only {len(stars)} stars detected")
            return None

        # Measure PSF for brightest stars
        fwhm_list = []
        for i, (x, y, bright) in enumerate(sorted(stars, key=lambda s: -s[2])[:num_stars]):
            fwhm = measure_psf_fwhm(image, x, y)
            if fwhm is not None and 0.5 < fwhm < 20:  # Sanity check
                fwhm_list.append(fwhm)

        if not fwhm_list:
            if verbose:
                print(f"  Warning: Could not fit any PSFs")
            return None

        # Average FWHM
        fwhm_pix = np.median(fwhm_list)
        pixel_scale = get_pixel_scale(instrument)
        seeing_arcsec = fwhm_pix * pixel_scale

        if verbose:
            print(f"  FWHM: {fwhm_pix:.2f} pix → {seeing_arcsec:.2f} arcsec")

        return seeing_arcsec

    except Exception as e:
        if verbose:
            print(f"  Error: {e}")
        return None


def measure_seeing_batch(filenames: Union[List[str], pd.Series],
                        instruments: Union[List[str], pd.Series],
                        threshold: float = 3.0,
                        verbose: bool = True,
                        use_header_keywords: Optional[Dict[str, str]] = None) -> np.ndarray:
    """
    Measure seeing for a batch of FITS files.

    For instruments where header keywords contain seeing measurements
    (e.g., LT/IOO with L1SEESEC), use those instead of measuring.

    Parameters
    ----------
    filenames : list or pd.Series
        List of FITS file paths
    instruments : list or pd.Series
        Corresponding instrument names
    threshold : float
        Star detection threshold (sigma)
    verbose : bool
        Print progress
    use_header_keywords : dict, optional
        Instrument -> header keyword mapping for pre-measured seeing.
        E.g., {'LT/IOO': 'L1SEESEC'}

    Returns
    -------
    seeing_values : ndarray
        Seeing in arcseconds for each file (NaN if measurement failed)
    """
    if use_header_keywords is None:
        use_header_keywords = {'LT/IOO': 'L1SEESEC'}

    seeing = np.full(len(filenames), np.nan)

    for i, (filename, instrument) in enumerate(zip(filenames, instruments)):
        if verbose:
            print(f"[{i+1}/{len(filenames)}] {Path(filename).name}...", end=' ')

        # Check if we should use header keyword instead of measuring
        if instrument in use_header_keywords:
            try:
                with fits.open(filename) as hdul:
                    keyword = use_header_keywords[instrument]
                    seeing[i] = float(hdul[0].header.get(keyword, np.nan))
                    if verbose:
                        if np.isnan(seeing[i]):
                            print(f"NO {keyword}")
                        else:
                            print(f"✓ {seeing[i]:.2f}\" ({keyword})\"")
            except Exception as e:
                if verbose:
                    print(f"ERROR: {e}")
        else:
            # Measure seeing from image
            seeing[i] = measure_seeing_file(filename, instrument, threshold=threshold,
                                           verbose=verbose)

            if verbose:
                if np.isnan(seeing[i]):
                    print("FAILED")
                else:
                    print(f"✓ {seeing[i]:.2f}\"")

    return seeing


def add_seeing_to_dataframe(df: pd.DataFrame,
                           threshold: float = 3.0,
                           verbose: bool = True) -> pd.DataFrame:
    """
    Add 'seeing' column to DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Must have 'filename' and 'instrument' columns
    threshold : float
        Star detection threshold
    verbose : bool
        Print progress

    Returns
    -------
    df_with_seeing : pd.DataFrame
        DataFrame with 'seeing' column added
    """
    required = ['filename', 'instrument']
    for col in required:
        if col not in df.columns:
            raise ValueError(f"DataFrame missing required column: '{col}'")

    df = df.copy()
    df['seeing'] = measure_seeing_batch(df['filename'], df['instrument'],
                                       threshold=threshold, verbose=verbose)

    return df


# Quick utility to summarize seeing stats
def seeing_summary(seeing: np.ndarray, by_filter: Optional[pd.Series] = None,
                   by_instrument: Optional[pd.Series] = None) -> dict:
    """
    Get summary statistics of seeing measurements.

    Parameters
    ----------
    seeing : ndarray
        Array of seeing values
    by_filter : pd.Series, optional
        Filter names for grouping
    by_instrument : pd.Series, optional
        Instrument names for grouping

    Returns
    -------
    summary : dict
        Statistics
    """
    valid = seeing[~np.isnan(seeing)]

    summary = {
        'count': len(valid),
        'mean': np.mean(valid),
        'median': np.median(valid),
        'std': np.std(valid),
        'min': np.min(valid),
        'max': np.max(valid),
        'failed': np.sum(np.isnan(seeing)),
    }

    if by_filter is not None:
        summary['by_filter'] = {}
        for filt in by_filter.unique():
            mask = by_filter == filt
            vals = seeing[mask]
            valid_vals = vals[~np.isnan(vals)]
            if len(valid_vals) > 0:
                summary['by_filter'][filt] = {
                    'median': np.median(valid_vals),
                    'mean': np.mean(valid_vals),
                    'count': len(valid_vals),
                }

    if by_instrument is not None:
        summary['by_instrument'] = {}
        for inst in by_instrument.unique():
            mask = by_instrument == inst
            vals = seeing[mask]
            valid_vals = vals[~np.isnan(vals)]
            if len(valid_vals) > 0:
                summary['by_instrument'][inst] = {
                    'median': np.median(valid_vals),
                    'mean': np.mean(valid_vals),
                    'count': len(valid_vals),
                }

    return summary
