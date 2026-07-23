"""
Color image creation from pandas DataFrame of multi-filter FITS files.

Designed for interactive notebook use. Handles WCS-based alignment and
intelligent RGB channel assignment.

Usage in notebook:
    from color_image_notebook import ColorImageMaker

    cim = ColorImageMaker(wavelengths={'g': 473, 'r': 622, 'i': 763})
    rgb, info = cim.create_from_dataframe(df)
    cim.plot(rgb)
"""

import warnings
from pathlib import Path
from typing import Dict, Tuple, Optional, List
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches

from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import ImageNormalize, AsinhStretch, PercentileInterval
from reproject import reproject_interp
from PIL import Image


class ColorImageMaker:
    """Create color composite images from multi-filter FITS files using WCS alignment."""

    def __init__(self, wavelengths: Optional[Dict[str, float]] = None,
                 scaling_method: str = 'asinh',
                 filter_keyword: str = 'FILTER',
                 verbose: bool = True):
        """
        Initialize the color image maker.

        Parameters
        ----------
        wavelengths : dict, optional
            Filter -> wavelength (nm) mapping. Defaults to standard filters.
        scaling_method : str
            'asinh', 'percentile', 'linear', or 'sqrt'
        filter_keyword : str
            FITS header keyword for filter name
        verbose : bool
            Print progress messages
        """
        self.verbose = verbose

        # Default wavelengths (SDSS/ZTF central wavelengths in nm)
        self.wavelengths = {
            'g': 473, 'r': 622, 'i': 763, 'z': 905,
            'y': 960, 'u': 355, 'B': 445, 'V': 551,
            'R': 658, 'I': 806
        }
        if wavelengths:
            self.wavelengths.update(wavelengths)

        self.scaling_method = scaling_method
        self.filter_keyword = filter_keyword
        self.images_raw = {}
        self.images_aligned = {}
        self.headers = {}
        self.wcs_solutions = {}
        self.rgb_mapping = {}
        self.rgb_image = None

    def _vprint(self, *args, **kwargs):
        """Print only if verbose."""
        if self.verbose:
            print(*args, **kwargs)

    def read_fits(self, filepath: str) -> Tuple[np.ndarray, fits.Header, WCS]:
        """
        Read FITS file and extract image, header, and WCS.

        Parameters
        ----------
        filepath : str
            Path to FITS file

        Returns
        -------
        image : ndarray
            2D image data
        header : fits.Header
            FITS header
        wcs : WCS
            WCS solution from header
        """
        with fits.open(filepath) as hdul:
            image = hdul[0].data
            header = hdul[0].header

        wcs = WCS(header)
        return image, header, wcs

    def extract_filter(self, header: fits.Header) -> str:
        """
        Extract filter name from FITS header.

        Parameters
        ----------
        header : fits.Header
            FITS header

        Returns
        -------
        filter_name : str
            Single character filter name
        """
        # Include FILTER1 for LT/IOO and other variants
        keywords = [self.filter_keyword, 'FILTER', 'FILTER1', 'FILT', 'FILTID', 'FILTERID', 'INSFLTR']

        for kw in keywords:
            if kw in header:
                filt = str(header[kw]).strip().lower()
                # Extract single character if embedded in longer string
                if len(filt) > 1:
                    # Check if it's in format like "BRG    " - take first char if it's a filter
                    if filt[0] in 'grizuybvri':
                        return filt[0]
                    # Otherwise look for single filter char
                    for char in filt:
                        if char in 'grizuybvri':
                            return char
                elif filt in 'grizuybvri':
                    return filt

        raise ValueError(f"Could not identify filter in header. Tried: {keywords}")

    def scale_image(self, image: np.ndarray, method: Optional[str] = None) -> np.ndarray:
        """
        Scale image to 0-255 range.

        Parameters
        ----------
        image : ndarray
            Raw image data
        method : str, optional
            Scaling method. Uses self.scaling_method if not specified.

        Returns
        -------
        scaled : ndarray
            Scaled image (0-255)
        """
        method = method or self.scaling_method
        valid_mask = np.isfinite(image)

        if not np.any(valid_mask):
            raise ValueError("Image contains no valid data")

        if method == 'asinh':
            norm = ImageNormalize(image, interval=PercentileInterval(99.5),
                                stretch=AsinhStretch(a=0.02))
        elif method == 'percentile':
            vmin, vmax = np.nanpercentile(image[valid_mask], [1, 99])
            norm = plt.Normalize(vmin=vmin, vmax=vmax, clip=True)
        elif method == 'sqrt':
            vmin, vmax = np.nanpercentile(image[valid_mask], [1, 99])
            image_sqrt = np.sqrt(np.clip(image, vmin, vmax) - vmin)
            return (image_sqrt / np.nanmax(image_sqrt) * 255).astype(np.uint8)
        elif method == 'linear':
            vmin, vmax = np.nanpercentile(image[valid_mask], [1, 99])
            norm = plt.Normalize(vmin=vmin, vmax=vmax, clip=True)
        else:
            raise ValueError(f"Unknown scaling method: {method}")

        scaled = norm(image)
        return (np.clip(scaled, 0, 1) * 255).astype(np.uint8)

    def align_to_reference(self, reference_wcs: WCS,
                          reference_shape: Tuple[int, int],
                          images_dict: Dict[str, np.ndarray],
                          wcs_dict: Dict[str, WCS]) -> Dict[str, np.ndarray]:
        """
        Align all images to a reference WCS using reproject_interp.

        Parameters
        ----------
        reference_wcs : WCS
            Target WCS solution
        reference_shape : tuple
            Target output shape (height, width)
        images_dict : dict
            Filter -> image array mapping
        wcs_dict : dict
            Filter -> WCS mapping

        Returns
        -------
        aligned : dict
            Filter -> aligned image array mapping
        """
        aligned = {}

        for filt, image in images_dict.items():
            self._vprint(f"  Aligning {filt}...", end=' ')

            if image.shape == reference_shape:
                # Check if WCS is already compatible
                if wcs_dict[filt].wcs.ctype == reference_wcs.wcs.ctype:
                    aligned[filt] = image
                    self._vprint("already aligned")
                    continue

            try:
                reprojected, footprint = reproject_interp(
                    (image, wcs_dict[filt]),
                    reference_wcs,
                    shape_out=reference_shape,
                    order='bicubic'
                )
                # Replace NaN pixels with 0
                reprojected = np.nan_to_num(reprojected, nan=0)
                aligned[filt] = reprojected
                self._vprint(f"✓ {reprojected.shape}")
            except Exception as e:
                self._vprint(f"✗ Error: {e}")
                raise

        return aligned

    def map_filters_to_rgb(self, filters: List[str]) -> Dict[str, str]:
        """
        Map filters to RGB based on wavelengths.

        Parameters
        ----------
        filters : list
            Available filter names

        Returns
        -------
        mapping : dict
            Filter -> 'R'/'G'/'B' mapping
        """
        if len(filters) < 3:
            raise ValueError(f"Need 3+ filters for RGB, got {len(filters)}: {filters}")

        # Sort by wavelength
        sorted_filters = sorted(filters, key=lambda f: self.wavelengths[f])

        mapping = {
            sorted_filters[-1]: 'R',  # Longest -> Red
            sorted_filters[1]: 'G',   # Middle -> Green
            sorted_filters[0]: 'B',   # Shortest -> Blue
        }

        self._vprint("\nFilter -> RGB mapping:")
        for filt in sorted(filters, key=lambda f: self.wavelengths[f]):
            channel = mapping[filt]
            print(f"  {filt} ({self.wavelengths[filt]:4.0f} nm) -> {channel}")

        return mapping

    def create_from_dataframe(self, df: pd.DataFrame,
                             output_path: Optional[str] = None,
                             check_mjd_range: Optional[float] = None) -> Tuple[np.ndarray, dict]:
        """
        Create color image from DataFrame of FITS files.

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame with columns: filter, mjd, instrument, filename
            Can have additional columns (ignored).
        output_path : str, optional
            Save output to this file (e.g., 'color.png')
        check_mjd_range : float, optional
            Check that all MJD values are within this range (days).
            If provided and violated, raise warning.

        Returns
        -------
        rgb_image : ndarray
            RGB image array (H x W x 3, uint8)
        info : dict
            Metadata about the image creation
        """
        # Validate DataFrame
        required_cols = ['filter', 'mjd', 'instrument', 'filename']
        for col in required_cols:
            if col not in df.columns:
                raise ValueError(f"DataFrame missing required column: '{col}'")

        self._vprint(f"Creating color image from {len(df)} files")
        print(f"Filters: {sorted(df['filter'].unique())}")
        print(f"MJD range: {df['mjd'].min():.3f} - {df['mjd'].max():.3f}")
        print(f"Instruments: {sorted(df['instrument'].unique())}")

        # Check MJD range
        if check_mjd_range:
            mjd_span = df['mjd'].max() - df['mjd'].min()
            if mjd_span > check_mjd_range:
                warnings.warn(
                    f"MJD span ({mjd_span:.3f} days) exceeds threshold "
                    f"({check_mjd_range} days). Images may not be well-aligned."
                )

        # Read files
        self._vprint("\nReading FITS files...")
        for idx, row in df.iterrows():
            filepath = row['filename']
            filt = row['filter'].lower()

            try:
                image, header, wcs = self.read_fits(filepath)
            except Exception as e:
                self._vprint(f"  ✗ {Path(filepath).name}: {e}")
                continue

            # Use filter from DataFrame (already extracted and validated)
            # Only warn if we can read it and it differs
            try:
                header_filt = self.extract_filter(header).lower()
                if filt != header_filt:
                    self._vprint(
                        f"  Warning: {Path(filepath).name} - "
                        f"DataFrame says {filt} but header says {header_filt}"
                    )
            except:
                # If we can't extract filter from header, that's OK - we'll use DataFrame value
                pass

            self.images_raw[filt] = image
            self.headers[filt] = header
            self.wcs_solutions[filt] = wcs
            self._vprint(f"  ✓ {Path(filepath).name}: {filt} {image.shape}")

        available_filters = list(self.images_raw.keys())
        self._vprint(f"\nLoaded filters: {available_filters}")

        # Get RGB mapping
        self.rgb_mapping = self.map_filters_to_rgb(available_filters)

        # Choose reference (first filter in wavelength order for stability)
        reference_filter = sorted(available_filters,
                                 key=lambda f: self.wavelengths[f])[1]  # Use middle filter
        reference_wcs = self.wcs_solutions[reference_filter]
        reference_shape = self.images_raw[reference_filter].shape

        self._vprint(f"\nUsing '{reference_filter}' as reference "
                     f"(shape: {reference_shape})")

        # Align images
        self._vprint("\nAligning images to reference WCS...")
        self.images_aligned = self.align_to_reference(
            reference_wcs, reference_shape,
            self.images_raw, self.wcs_solutions
        )

        # Create RGB
        self._vprint("\nScaling and combining into RGB...")
        rgb_image = np.zeros((*reference_shape, 3), dtype=np.uint8)
        channel_map = {'R': 0, 'G': 1, 'B': 2}

        for filt in available_filters:
            channel = self.rgb_mapping[filt]
            channel_idx = channel_map[channel]
            scaled = self.scale_image(self.images_aligned[filt])
            rgb_image[:, :, channel_idx] = scaled
            self._vprint(f"  {filt} -> {channel}: range [{scaled.min()}, {scaled.max()}]")

        self.rgb_image = rgb_image

        # Save if requested
        if output_path:
            img = Image.fromarray(rgb_image, mode='RGB')
            img.save(output_path)
            self._vprint(f"\n✓ Saved to {output_path}")

        info = {
            'filters': available_filters,
            'rgb_mapping': self.rgb_mapping,
            'reference_filter': reference_filter,
            'shape': reference_shape,
            'scaling_method': self.scaling_method,
            'num_images': len(df),
        }

        return rgb_image, info

    def plot(self, rgb_image: Optional[np.ndarray] = None,
             figsize: Tuple[int, int] = (12, 10),
             title: str = 'Color Composite Image',
             show_info: bool = True) -> plt.Figure:
        """
        Display the color image in a notebook.

        Parameters
        ----------
        rgb_image : ndarray, optional
            Image to plot. Uses self.rgb_image if not provided.
        figsize : tuple
            Figure size (width, height)
        title : str
            Plot title
        show_info : bool
            Show filter->channel mapping info

        Returns
        -------
        fig : matplotlib.figure.Figure
            Matplotlib figure object
        """
        rgb = rgb_image if rgb_image is not None else self.rgb_image

        if rgb is None:
            raise ValueError("No image to plot. Run create_from_dataframe() first.")

        fig, ax = plt.subplots(figsize=figsize)
        ax.imshow(rgb)
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.axis('off')

        # Add info text
        if show_info and self.rgb_mapping:
            info_text = "Channels:\n"
            for filt, channel in sorted(self.rgb_mapping.items(),
                                        key=lambda x: ['B', 'G', 'R'].index(x[1])):
                wl = self.wavelengths.get(filt, '?')
                info_text += f"  {channel}: {filt} ({wl} nm)\n"

            ax.text(0.02, 0.98, info_text, transform=ax.transAxes,
                   fontsize=10, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
                   family='monospace')

        plt.tight_layout()
        return fig

    def save(self, output_path: str) -> None:
        """
        Save the current RGB image to file.

        Parameters
        ----------
        output_path : str
            Output filename (e.g., 'color.png', 'color.jpg')
        """
        if self.rgb_image is None:
            raise ValueError("No image to save. Run create_from_dataframe() first.")

        img = Image.fromarray(self.rgb_image, mode='RGB')
        img.save(output_path)
        self._vprint(f"✓ Saved to {output_path}")

    def get_info(self) -> dict:
        """Get metadata about the current image."""
        if not self.rgb_mapping:
            return {}

        return {
            'rgb_mapping': self.rgb_mapping,
            'wavelengths': self.wavelengths,
            'scaling_method': self.scaling_method,
            'shape': self.rgb_image.shape if self.rgb_image is not None else None,
        }
