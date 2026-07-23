"""
Colour composite maker for multi-filter / multi-telescope astronomical imaging.

Pipeline (all on a common WCS centred on the target):

    1. Reproject every input frame onto a TAN WCS centred on ``center_coord``
       using the Drizzle-style adaptive reprojection (no moiré).
    2. Per-frame: sigma-clipped sky removal and EXPTIME normalisation.
    3. Stack same-filter frames (sigma-clipped exposure-time-weighted mean).
    4. Source-masked Gaussian large-scale background subtraction.
    5. Star detection + cross-matched sub-pixel shift to align the three
       filters to the reference (no colour halos around bright stars).
    6. PSF homogenisation: convolve sharper filters to match the broadest
       per-filter FWHM.
    7. Per-channel scaling to a common bright-pixel level, optional weighting,
       then combine with ``astropy.visualization.make_lupton_rgb``.

Public entry points:
    make_color_composite(df, ...)         -- one RGB composite from a DataFrame
    find_nightly_groups(df, ...)          -- group rows by night (MJD - 0.5)
    find_temporal_groups(df, ...)         -- group rows by a sliding ±N-day window
    save_color_pdf(rgb, info, ...)        -- PDF with target circle + zoom
    make_nightly_pdfs(df, ...)            -- one PDF per night
    make_temporal_group_pdfs(df, ...)     -- one PDF per temporal window
    make_color_from_dataframe(df, ...)    -- batch helper used by older notebooks
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import SigmaClip, sigma_clipped_stats, sigma_clip
from astropy.visualization import make_lupton_rgb
from astropy.wcs import WCS, FITSFixedWarning
from PIL import Image
from photutils.background import Background2D, MedianBackground, SExtractorBackground
from photutils.detection import DAOStarFinder
from reproject import reproject_interp, reproject_adaptive
from scipy.ndimage import gaussian_filter, shift as ndi_shift

warnings.simplefilter("ignore", FITSFixedWarning)

# Optional cosmic-ray rejection (LACosmic via astroscrappy). Kept optional so
# the module still imports cleanly on machines without astroscrappy installed.
try:
    from astroscrappy import detect_cosmics as _detect_cosmics
    _HAS_ASTROSCRAPPY = True
except ImportError:  # pragma: no cover
    _detect_cosmics = None
    _HAS_ASTROSCRAPPY = False


DEFAULT_WAVELENGTHS = {
    "u": 355,
    "g": 473,
    "r": 622,
    "i": 763,
    "z": 905,
    "y": 960,
}


# ---------------------------------------------------------------------------
# Logging helper
# ---------------------------------------------------------------------------

def _vprint(verbose: bool, *args, **kwargs) -> None:
    if verbose:
        print(*args, **kwargs)


# ---------------------------------------------------------------------------
# FITS / WCS helpers
# ---------------------------------------------------------------------------

def _open_image(path: str) -> Tuple[np.ndarray, WCS, fits.Header]:
    """Open the first HDU that has 2-D image data and a usable celestial WCS."""
    with fits.open(path, memmap=False) as hdul:
        for hdu in hdul:
            if hdu.data is None or hdu.data.ndim != 2:
                continue
            try:
                wcs = WCS(hdu.header).celestial
                if wcs.has_celestial:
                    return np.array(hdu.data, dtype=float), wcs, hdu.header
            except Exception:
                continue
    raise ValueError(f"No 2-D image with WCS in {path}")


def _pixel_scale_arcsec(wcs: WCS) -> float:
    scales = wcs.proj_plane_pixel_scales()
    return float(np.mean([s.to("arcsec").value for s in scales]))


# ---------------------------------------------------------------------------
# Cosmic-ray rejection (LACosmic / astroscrappy)
# ---------------------------------------------------------------------------

# Sensible default detect_cosmics parameters. These work well for typical
# ground-based optical CCD frames (NOT/ALFOSC, LT/IOO, TJO/MEIA2 etc.).
# Overridable per call via the ``lacosmic_kwargs`` argument of
# ``make_color_composite``.
_LACOSMIC_DEFAULTS = dict(
    sigclip=4.5,        # Laplacian S/N threshold for a pixel to be flagged
    sigfrac=0.3,        # neighbour-pixel growth threshold (fraction of sigclip)
    objlim=5.0,         # contrast vs. fine-structure image (raise to spare PSFs)
    readnoise=5.0,      # e- RMS read noise (overridden by header RDNOISE if present)
    gain=1.0,           # e-/ADU (overridden by header GAIN if present)
    satlevel=60000.0,   # ADU saturation level (overridden by header SATURATE)
    niter=4,
    cleantype="meanmask",
    fsmode="median",
    verbose=False,
)


def _clean_cosmics(data: np.ndarray,
                   header: fits.Header,
                   lacosmic_kwargs: Optional[Dict] = None,
                   ) -> Tuple[np.ndarray, int]:
    """Run LACosmic (astroscrappy) on a single frame.

    Parameters
    ----------
    data
        2-D float array of the science frame, in ADU.
    header
        FITS header — used to look up GAIN, RDNOISE and SATURATE so the
        detector-specific values override the generic defaults.
    lacosmic_kwargs
        Optional overrides merged on top of ``_LACOSMIC_DEFAULTS``.

    Returns
    -------
    cleaned : 2-D array with CRs replaced by the local median (same shape, same dtype).
    n_cr    : number of pixels flagged as cosmic rays.

    If astroscrappy is not installed, the input array is returned unchanged
    with ``n_cr = -1`` so the caller can warn appropriately.
    """
    if not _HAS_ASTROSCRAPPY:
        return data, -1

    params = dict(_LACOSMIC_DEFAULTS)
    # Pull detector-specific values from the header when available so the
    # CR detection threshold scales with the noise floor of each instrument.
    for hdr_key, kw_key in (("GAIN", "gain"),
                            ("EGAIN", "gain"),
                            ("RDNOISE", "readnoise"),
                            ("READNOIS", "readnoise"),
                            ("RON", "readnoise"),
                            ("SATURATE", "satlevel"),
                            ("SATLEVEL", "satlevel")):
        if hdr_key in header:
            try:
                val = float(header[hdr_key])
                if np.isfinite(val) and val > 0:
                    params[kw_key] = val
            except (TypeError, ValueError):
                pass
    if lacosmic_kwargs:
        params.update(lacosmic_kwargs)

    # astroscrappy needs a finite float32 array. NaNs in the input would
    # otherwise propagate into the Laplacian and the result.
    work = np.ascontiguousarray(np.nan_to_num(data, nan=0.0), dtype=np.float32)

    try:
        mask, cleaned = _detect_cosmics(work, **params)
    except Exception as e:
        # If astroscrappy fails (e.g. extreme params, all-NaN frame), fall
        # back to the original data rather than crashing the whole pipeline.
        warnings.warn(f"astroscrappy.detect_cosmics failed: {e}; "
                      "returning input frame unchanged.", RuntimeWarning)
        return data, 0

    # Preserve NaN locations from the input (we filled them with 0 for
    # astroscrappy) so downstream sky estimation still ignores them.
    cleaned = cleaned.astype(data.dtype, copy=False)
    nan_in = ~np.isfinite(data)
    if nan_in.any():
        cleaned = np.where(nan_in, np.nan, cleaned)

    return cleaned, int(mask.sum())


def _build_target_wcs(center: SkyCoord, pixscale_arcsec: float, size: int) -> WCS:
    w = WCS(naxis=2)
    w.wcs.crpix = [(size + 1) / 2.0, (size + 1) / 2.0]
    w.wcs.crval = [center.ra.deg, center.dec.deg]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.cdelt = [-pixscale_arcsec / 3600.0, pixscale_arcsec / 3600.0]
    w.wcs.cunit = ["deg", "deg"]
    return w


def _normalise_filter(name: str) -> str:
    s = str(name).strip().lower()
    for c in s:
        if c in "ugrizy":
            return c
    return s[:1]


def _sigma_clipped_sky(image: np.ndarray) -> Tuple[float, float]:
    finite = np.isfinite(image)
    if not finite.any():
        return 0.0, 0.0
    _, median, std = sigma_clipped_stats(image[finite], sigma=3.0, maxiters=5)
    return float(median), float(std)


# ---------------------------------------------------------------------------
# Reproject + per-frame sky subtraction
# ---------------------------------------------------------------------------

def _load_and_reproject(paths: List[str],
                        target_wcs: WCS,
                        target_shape: Tuple[int, int],
                        verbose: bool,
                        lacosmic: bool = False,
                        lacosmic_kwargs: Optional[Dict] = None,
                        ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """
    Reproject every frame onto the target grid.

    Parameters
    ----------
    lacosmic : bool, optional (default False)
        If True, run astroscrappy.detect_cosmics on each raw frame *before*
        reprojection. Cosmic-ray-flagged pixels are replaced by a local
        median so the subsequent reproject + stack produces a clean image.
        Requires the optional ``astroscrappy`` package; if unavailable the
        step is skipped with a warning.
    lacosmic_kwargs : dict, optional
        Per-call overrides for the detect_cosmics parameters
        (e.g. ``{"sigclip": 5.0, "objlim": 6.0}``). Merged on top of
        ``_LACOSMIC_DEFAULTS`` and any GAIN / RDNOISE / SATURATE values
        read from each frame's FITS header.

    Returns
    -------
    layers : (n_frames, H, W) float array  (NaN outside the source frame)
    weights : (n_frames,) float array  (per-frame exposure time, default 1)
    coverage : (H, W) integer array
    total_exptime : float
    """
    layers = []
    weights = []
    coverage = np.zeros(target_shape, dtype=np.int32)
    total_exp = 0.0

    if lacosmic and not _HAS_ASTROSCRAPPY:
        warnings.warn("lacosmic=True requested but astroscrappy is not "
                      "installed; cosmic-ray rejection step skipped. "
                      "Install with: pip install astroscrappy",
                      RuntimeWarning)
        lacosmic = False

    for path in paths:
        try:
            data, wcs, header = _open_image(path)
        except Exception as e:
            _vprint(verbose, f"      skip {Path(path).name}: {e}")
            continue

        # ----------------------------------------------------------------
        # Cosmic-ray rejection (per raw frame, before reprojection).
        # Doing this pre-reprojection keeps the CR signatures sharp and
        # easy for LACosmic's Laplacian edge detector to find. Once
        # frames are resampled onto the target grid the CR cores are
        # spread over neighbours and become much harder to clean.
        # ----------------------------------------------------------------
        if lacosmic:
            data, n_cr = _clean_cosmics(data, header,
                                        lacosmic_kwargs=lacosmic_kwargs)
            if n_cr > 0:
                _vprint(verbose, f"      lacosmic: cleaned {n_cr:6d} "
                                 f"CR pixels in {Path(path).name}")

        exp = float(header.get("EXPTIME", 0.0) or 0.0)

        try:
            # reproject_adaptive (Drizzle-style) handles oversampling without
            # the moire / grid pattern that bilinear interp produces when going
            # from a coarser pixel grid (e.g. LT 0.30"/px) onto a finer one
            # (e.g. NOT 0.21"/px).
            arr, fp = reproject_adaptive(
                (np.nan_to_num(data, nan=0.0), wcs),
                target_wcs,
                shape_out=target_shape,
                kernel="gaussian",
                conserve_flux=False,
            )
        except Exception as e:
            _vprint(verbose, f"      adaptive reproject failed on "
                              f"{Path(path).name} ({e}), trying bilinear...")
            try:
                arr, fp = reproject_interp(
                    (np.nan_to_num(data, nan=0.0), wcs),
                    target_wcs,
                    shape_out=target_shape,
                    order="bilinear",
                )
            except Exception as e2:
                _vprint(verbose, f"      reproject failed on {Path(path).name}: {e2}")
                continue

        mask = fp.astype(bool)
        if mask.sum() == 0:
            continue

        # Per-frame sky subtraction so frames combine on a common pedestal
        sky, _ = _sigma_clipped_sky(arr[mask])
        arr = arr - sky
        if exp > 0:
            arr = arr / exp
            total_exp += exp
            w_frame = exp
        else:
            w_frame = 1.0

        layers.append(np.where(mask, arr, np.nan))
        weights.append(w_frame)
        coverage += mask.astype(np.int32)

    if not layers:
        raise RuntimeError("No frames could be reprojected")

    return (np.stack(layers, axis=0),
            np.asarray(weights, dtype=float),
            coverage,
            total_exp)


# ---------------------------------------------------------------------------
# Stacking: sigma-clipped exposure-time-weighted mean
# ---------------------------------------------------------------------------

def _stack_weighted(layers: np.ndarray,
                    weights: np.ndarray,
                    sigma: float = 3.0) -> np.ndarray:
    """
    Sigma-clipped weighted mean stack along axis 0.

    NaN pixels are ignored. Pixels rejected by sigma_clip are also ignored.
    """
    # Use sigma_clip *across* the stack axis so cosmic rays / outliers in any
    # single frame are masked.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        masked = sigma_clip(layers, sigma=sigma, axis=0, masked=True, copy=False)

    # masked is a MaskedArray with shape == layers.shape
    data = masked.data.astype(float)
    valid = (~masked.mask) & np.isfinite(data)

    w = np.asarray(weights, dtype=float).reshape(-1, 1, 1)
    w = np.broadcast_to(w, data.shape) * valid

    sum_w = w.sum(axis=0)
    sum_wx = (data * w).sum(axis=0)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        out = np.where(sum_w > 0, sum_wx / sum_w, 0.0)
    return out


# ---------------------------------------------------------------------------
# 2-D background subtraction
# ---------------------------------------------------------------------------

def _subtract_background_2d(image: np.ndarray,
                            large_scale_sigma_px: float = 60.0,
                            box_size: int = 64,           # accepted for API compatibility
                            filter_size: int = 5) -> Tuple[np.ndarray, float]:
    """Fit + subtract a slowly varying sky map; return (residual, background_median).

    Pure Gaussian-smoothed background model.

    Earlier versions used ``photutils.Background2D`` plus a Gaussian
    large-scale component, but the box-grid produced subtle 64-px hue
    offsets on the final colour image because each filter had slightly
    different per-box sky estimates. A single source-masked Gaussian
    smoothing is grid-free and gives equally good slow-component removal
    in practice for this use case.

    Procedure
    ---------
      1. Sigma-clipped global sky (gives the noise floor and a fallback sky).
      2. Source mask: pixels ≥5σ above the noise floor, dilated so PSF wings
         and bright extended emission are excluded but diffuse haze isn't.
      3. Replace masked pixels with the global sky.
      4. Gaussian smooth (σ ≈ 60 px) → background model.
      5. Subtract.
    """
    finite = np.isfinite(image)
    if not finite.any():
        return image, 0.0

    _, sky0, sky_std0 = sigma_clipped_stats(image[finite], sigma=3.0)
    if sky_std0 <= 0:
        sky_std0 = float(np.nanstd(image[finite]))

    raw_mask = (image - sky0) > 5.0 * sky_std0
    src_mask = gaussian_filter(raw_mask.astype(float), sigma=3.0) > 0.05
    src_mask |= ~finite

    filled = np.where(src_mask, sky0, image)
    bg = gaussian_filter(filled, sigma=large_scale_sigma_px)

    return image - bg, float(np.nanmedian(bg))


# ---------------------------------------------------------------------------
# Star detection + sub-pixel registration refinement
# ---------------------------------------------------------------------------

def _detect_stars(image: np.ndarray, fwhm: float = 4.0,
                  threshold_sigma: float = 8.0,
                  max_stars: int = 80,
                  min_stars: int = 5,
                  fallback_thresholds=(5.0, 3.5, 2.5)) -> np.ndarray:
    """
    Run DAOStarFinder on a sky-subtracted image. Returns (N, 2) array of
    (x, y) centroids of the brightest non-saturated detections.

    If the requested ``threshold_sigma`` finds fewer than ``min_stars``,
    automatically retry at the lower thresholds in ``fallback_thresholds``.
    This is essential for deep stacked reference images (PS1, Legacy Survey)
    where 8σ is often too strict.
    """
    finite = np.isfinite(image)
    if not finite.any():
        return np.empty((0, 2))
    _, _, std = sigma_clipped_stats(image[finite], sigma=3.0)
    if std <= 0:
        return np.empty((0, 2))

    image_safe = np.where(finite, image, 0.0)
    thresholds = [float(threshold_sigma)] + [float(t) for t in fallback_thresholds]
    best = np.empty((0, 2))
    for t in thresholds:
        try:
            finder = DAOStarFinder(fwhm=fwhm, threshold=t * std,
                                   exclude_border=True)
            tbl = finder(image_safe)
        except Exception:
            tbl = None
        if tbl is None or len(tbl) == 0:
            continue
        tbl.sort("flux", reverse=True)
        tbl = tbl[:max_stars]
        coords = np.column_stack([np.asarray(tbl["xcentroid"]),
                                  np.asarray(tbl["ycentroid"])])
        best = coords
        if len(coords) >= min_stars:
            break
    return best


def _measure_fwhm(image: np.ndarray, fwhm_init: float = 4.0,
                  threshold_sigma: float = 10.0) -> float:
    """
    Estimate the median FWHM in pixels from the brightest few stars by fitting
    a 2-D Gaussian-equivalent quadratic to small cutouts. Robust enough for
    visualisation / convolution kernel sizing.

    Uses ``_detect_stars`` adaptive fallback so this works on both bright
    science stacks and deeper / lower-SNR reference stacks.
    """
    coords = _detect_stars(image, fwhm=fwhm_init,
                           threshold_sigma=threshold_sigma, max_stars=40)
    if len(coords) < 5:
        return float(fwhm_init)
    half = max(int(round(fwhm_init * 3)), 6)
    fwhms = []
    h, w = image.shape
    for x, y in coords:
        ix, iy = int(round(x)), int(round(y))
        if ix - half < 0 or ix + half >= w or iy - half < 0 or iy + half >= h:
            continue
        cut = image[iy - half:iy + half + 1, ix - half:ix + half + 1]
        if not np.all(np.isfinite(cut)):
            continue
        cut = cut - np.nanmedian(cut)
        peak = cut.max()
        if peak <= 0:
            continue
        # Half-maximum width along x and y (linear interpolation)
        cx, cy = cut.shape[1] // 2, cut.shape[0] // 2
        try:
            xprof = cut[cy, :]
            yprof = cut[:, cx]
            def _hw(prof, peak):
                above = prof >= peak / 2.0
                if not above.any():
                    return None
                idx = np.where(above)[0]
                return idx.max() - idx.min() + 1
            hx = _hw(xprof, peak)
            hy = _hw(yprof, peak)
            if hx and hy:
                fwhms.append(0.5 * (hx + hy))
        except Exception:
            continue
    if not fwhms:
        return float(fwhm_init)
    return float(np.median(fwhms))


def _match_offset(stars_ref: np.ndarray,
                  stars_other: np.ndarray,
                  search_radius_px: float = 6.0) -> Tuple[float, float, int]:
    """
    Estimate the (dx, dy) pixel shift that brings ``stars_other`` onto
    ``stars_ref``. Robust to wrong matches via median over nearest-neighbour
    differences with sigma clipping.
    """
    if len(stars_ref) < 3 or len(stars_other) < 3:
        return 0.0, 0.0, 0

    diffs = []
    for x, y in stars_other:
        d = np.hypot(stars_ref[:, 0] - x, stars_ref[:, 1] - y)
        j = int(np.argmin(d))
        if d[j] <= search_radius_px:
            diffs.append((stars_ref[j, 0] - x, stars_ref[j, 1] - y))
    if len(diffs) < 3:
        return 0.0, 0.0, 0
    diffs = np.asarray(diffs)
    dx = sigma_clipped_stats(diffs[:, 0], sigma=2.5)[1]
    dy = sigma_clipped_stats(diffs[:, 1], sigma=2.5)[1]
    return float(dx), float(dy), len(diffs)


# ---------------------------------------------------------------------------
# Main API
# ---------------------------------------------------------------------------

def _wavelength_to_rgb_weights(filter_wavelengths: Dict[str, float]) -> Dict[str, Tuple[float, float, float]]:
    """
    Triangular interpolation: each filter gets a (R, G, B) weight tuple
    based on its position in the wavelength range spanned by all filters.

    Endpoints get pure red / pure blue, the midpoint goes to pure green,
    intermediate filters split between two channels in proportion to their
    distance to the centre.
    """
    if not filter_wavelengths:
        return {}
    wls = list(filter_wavelengths.values())
    wl_min, wl_max = min(wls), max(wls)
    if wl_min >= wl_max:
        return {f: (1 / 3, 1 / 3, 1 / 3) for f in filter_wavelengths}
    wl_mid = 0.5 * (wl_min + wl_max)

    weights = {}
    for f, w in filter_wavelengths.items():
        if w <= wl_mid:
            t = (w - wl_min) / (wl_mid - wl_min)
            r, g, b = 0.0, t, 1.0 - t
        else:
            t = (w - wl_mid) / (wl_max - wl_mid)
            r, g, b = t, 1.0 - t, 0.0
        weights[f] = (r, g, b)
    return weights


def make_color_composite(df: pd.DataFrame,
                         filter_bands=None,                      # list[str] | "all" | None
                         output_path: Optional[str] = None,
                         wavelengths: Optional[Dict[str, float]] = None,
                         center_coord: Optional[str] = None,
                         cutout_size: int = 1024,
                         pixscale_arcsec: Optional[float] = None,
                         stretch: Optional[float] = None,
                         Q: float = 8.0,
                         minimum: float = 0.0,
                         channel_weights: Optional[Dict[str, float]] = None,
                         bkg_box_size: int = 64,
                         bkg_filter_size: int = 5,
                         match_psf: bool = True,
                         refine_alignment: bool = True,
                         lacosmic: bool = True,
                         lacosmic_kwargs: Optional[Dict] = None,
                         verbose: bool = True) -> Tuple[np.ndarray, Dict]:
    """
    Build a colour composite from a DataFrame of multi-filter FITS frames.

    Parameters
    ----------
    filter_bands : list[str] | "all" | None
        - ``None`` (default): the three reddest filters available, mapped
          R / G / B in order of decreasing wavelength.
        - ``list[str]``: three explicit filters in (R, G, B) order.
        - ``"all"``: use **every** filter present in the DataFrame, blending
          them into the R / G / B output channels by triangular wavelength
          weights (every filter contributes; longest-wavelength filters get
          pure R, shortest get pure B, intermediate filters split between
          two adjacent channels).
    lacosmic : bool, default True
        Run astroscrappy's LACosmic on every raw frame before reprojection
        so cosmic rays don't survive the stack as single-channel speckles
        in the final RGB composite. Requires the ``astroscrappy`` package;
        if not installed, the step is skipped with a warning. Set ``False``
        on archival reference cutouts (PS1, Legacy Survey) — they're
        already coadds and CR rejection there can clip real PSF cores.
    lacosmic_kwargs : dict, optional
        Per-call overrides for detect_cosmics parameters, e.g.
        ``{"sigclip": 5.0, "objlim": 6.0}``. The defaults
        (``sigclip=4.5, objlim=5.0``) are conservative enough not to
        clip stellar PSFs in well-sampled frames; raise ``objlim`` if
        you see bright stars losing flux.

    All the heavy lifting (alignment, stacking, sky removal, PSF matching,
    Lupton combine) is described in the module docstring.
    """
    if wavelengths is None:
        wavelengths = DEFAULT_WAVELENGTHS
    if center_coord is None:
        raise ValueError("center_coord is required (e.g. '14:37:16.14 +71:50:30.31')")

    # ------------------------------------------------------------------
    # 1. Group frames by filter
    # ------------------------------------------------------------------
    df = df.copy()
    df["_filt"] = df["filter"].apply(_normalise_filter)
    by_filter: Dict[str, List[str]] = {
        f: list(g["filename"]) for f, g in df.groupby("_filt")
    }
    available = sorted(by_filter.keys(), key=lambda f: wavelengths.get(f, 0))

    _vprint(verbose, f"Creating colour composite from {len(df)} frames")
    _vprint(verbose, f"Filters available: {available}  "
                     f"(counts: {[len(by_filter[f]) for f in available]})")
    if len(available) < 3:
        raise ValueError(f"Need at least 3 filters, got {available}")

    # Determine processing mode
    multifilter = (isinstance(filter_bands, str) and filter_bands.lower() == "all")
    if multifilter:
        active_filters = list(available)                # use everything
        rgb_weights = _wavelength_to_rgb_weights(
            {f: float(wavelengths.get(f, 0)) for f in active_filters})
        _vprint(verbose, "Multi-filter blend mode: every filter contributes.")
        for f in active_filters:
            r, g, b = rgb_weights[f]
            _vprint(verbose, f"  {f} ({wavelengths.get(f, '?')} nm) → "
                              f"R={r:.2f} G={g:.2f} B={b:.2f}")
    else:
        if filter_bands is None:
            filter_bands = [available[-1], available[-2], available[-3]]
        else:
            filter_bands = list(filter_bands)
        for f in filter_bands:
            if f not in by_filter:
                raise ValueError(f"Filter '{f}' not present in the DataFrame")
        active_filters = list(filter_bands)
        # Pure R/G/B for the three named filters
        rgb_weights = {
            filter_bands[0]: (1.0, 0.0, 0.0),
            filter_bands[1]: (0.0, 1.0, 0.0),
            filter_bands[2]: (0.0, 0.0, 1.0),
        }
        _vprint(verbose, f"Filter mapping  R={filter_bands[0]}  "
                         f"G={filter_bands[1]}  B={filter_bands[2]}")

    # ------------------------------------------------------------------
    # 2. Build the reference (target) WCS
    # ------------------------------------------------------------------
    center = SkyCoord(center_coord, unit=("hourangle", "deg"))
    if pixscale_arcsec is None:
        scales = []
        for f in active_filters:
            try:
                _, w, _ = _open_image(by_filter[f][0])
                scales.append(_pixel_scale_arcsec(w))
            except Exception:
                pass
        pixscale_arcsec = float(min(scales)) if scales else 0.3
    target_wcs = _build_target_wcs(center, pixscale_arcsec, cutout_size)
    target_shape = (cutout_size, cutout_size)
    fov_arcmin = cutout_size * pixscale_arcsec / 60.0
    _vprint(verbose, f"Target WCS: TAN, pixscale={pixscale_arcsec:.3f}\"/px, "
                     f"size={cutout_size}×{cutout_size}px (FOV ≈ {fov_arcmin:.2f}')")
    _vprint(verbose, f"Centre: RA={center.ra.deg:.5f}°, Dec={center.dec.deg:.5f}°\n")

    # ------------------------------------------------------------------
    # 3. Reproject + sigma-clipped weighted mean stack per filter
    # ------------------------------------------------------------------
    stacked: Dict[str, np.ndarray] = {}
    n_frames: Dict[str, int] = {}
    exptimes: Dict[str, float] = {}
    for f in active_filters:
        _vprint(verbose, f"  [{f}] reprojecting {len(by_filter[f])} frames...")
        layers, w, _, total_exp = _load_and_reproject(
            by_filter[f], target_wcs, target_shape, verbose=verbose,
            lacosmic=lacosmic, lacosmic_kwargs=lacosmic_kwargs)
        stacked[f] = _stack_weighted(layers, w, sigma=3.0)
        n_frames[f] = layers.shape[0]
        exptimes[f] = total_exp
        _vprint(verbose, f"  [{f}] σ-clipped weighted-mean stack of "
                         f"{layers.shape[0]} frames (Σ exptime={total_exp:.1f}s)")

    # ------------------------------------------------------------------
    # 4. 2-D background subtraction
    # ------------------------------------------------------------------
    _vprint(verbose, "\nFitting + subtracting 2-D background per filter...")
    bg_levels: Dict[str, float] = {}
    for f in active_filters:
        stacked[f], bg_med = _subtract_background_2d(
            stacked[f], box_size=bkg_box_size, filter_size=bkg_filter_size)
        bg_levels[f] = bg_med
        _vprint(verbose, f"  [{f}] background median ≈ {bg_med:.4g}")

    # ------------------------------------------------------------------
    # 5. PSF measurement + sub-pixel WCS refinement
    # ------------------------------------------------------------------
    fwhms = {f: _measure_fwhm(stacked[f], fwhm_init=4.0) for f in active_filters}
    _vprint(verbose, "\nMeasured FWHM (px):  " +
            "  ".join(f"{f}={fwhms[f]:.2f}" for f in active_filters))

    shifts: Dict[str, Tuple[float, float, int]] = {f: (0.0, 0.0, 0) for f in active_filters}
    if refine_alignment:
        # Detect stars in each filter, pick whichever has the MOST stars as
        # the alignment reference. Falls back to the reddest filter if no
        # filter has enough detections.
        per_filter_stars = {}
        for f in active_filters:
            per_filter_stars[f] = _detect_stars(stacked[f], fwhm=max(3.0, fwhms[f]))
        # Pick filter with the most stars (ties broken by wavelength → reddest)
        ref = max(active_filters,
                  key=lambda f: (len(per_filter_stars[f]),
                                 wavelengths.get(f, 0)))
        ref_stars = per_filter_stars[ref]
        _vprint(verbose, f"\nRefining alignment to '{ref}' "
                         f"({len(ref_stars)} reference stars)...")
        for f in active_filters:
            if f == ref:
                continue
            other_stars = per_filter_stars[f]
            dx, dy, n = _match_offset(ref_stars, other_stars,
                                      search_radius_px=max(8.0, 1.5 * max(fwhms.values())))
            shifts[f] = (dx, dy, n)
            if n >= 3 and (abs(dx) > 0.05 or abs(dy) > 0.05):
                stacked[f] = ndi_shift(np.nan_to_num(stacked[f], nan=0.0),
                                       shift=(dy, dx), order=3, mode="nearest")
            _vprint(verbose, f"  [{f}] shift=({dx:+.2f}, {dy:+.2f}) px  ({n} matches)")

    # ------------------------------------------------------------------
    # 6. PSF homogenisation: convolve sharper filters to the broadest FWHM
    # ------------------------------------------------------------------
    if match_psf:
        target_fwhm = max(fwhms.values())
        _vprint(verbose, f"\nMatching PSFs to FWHM = {target_fwhm:.2f} px...")
        for f in active_filters:
            sigma_native = fwhms[f] / 2.3548
            sigma_target = target_fwhm / 2.3548
            extra = sigma_target ** 2 - sigma_native ** 2
            if extra > 0.01:
                sigma_kernel = np.sqrt(extra)
                stacked[f] = gaussian_filter(stacked[f], sigma=sigma_kernel)
                _vprint(verbose, f"  [{f}] convolved with σ={sigma_kernel:.2f} px")
            else:
                _vprint(verbose, f"  [{f}] already at target PSF")

    # ------------------------------------------------------------------
    # 7. Final sky removal + per-filter scaling to a common bright-pixel
    #    level (so colours of point sources are preserved).
    # ------------------------------------------------------------------
    weights = channel_weights or {f: 1.0 for f in active_filters}
    sky_stats: Dict[str, Tuple[float, float]] = {}
    bright_levels: Dict[str, float] = {}

    h_, _ = stacked[active_filters[0]].shape
    s0, s1 = int(0.2 * h_), int(0.8 * h_)

    final = {}
    for f in active_filters:
        img = stacked[f]
        sky, std = _sigma_clipped_sky(img)
        img = img - sky

        sample = img[s0:s1, s0:s1]
        sample = sample[np.isfinite(sample)]
        if sample.size:
            level = float(np.percentile(sample, 99.7))
        else:
            level = 1.0
        if level <= 0:
            level = max(std, 1e-6)
        bright_levels[f] = level

        img = img / level
        img = img * float(weights.get(f, 1.0))
        sky_stats[f] = (sky, std)
        final[f] = img

    # ------------------------------------------------------------------
    # 8. Combine filters into R / G / B using rgb_weights
    #
    #    For 3-filter mode each weight is (1, 0, 0) / (0, 1, 0) / (0, 0, 1).
    #    For 'all' mode, every filter contributes proportional to its
    #    triangular wavelength weight. Per-channel total weight is computed
    #    so heavily-loaded channels don't end up brighter than sparse ones.
    # ------------------------------------------------------------------
    H, W = final[active_filters[0]].shape
    R = np.zeros((H, W), dtype=float)
    G = np.zeros((H, W), dtype=float)
    B = np.zeros((H, W), dtype=float)
    wR = wG = wB = 0.0
    for f in active_filters:
        rw, gw, bw = rgb_weights[f]
        if rw: R += rw * final[f]; wR += rw
        if gw: G += gw * final[f]; wG += gw
        if bw: B += bw * final[f]; wB += bw
    if wR > 0: R /= wR
    if wG > 0: G /= wG
    if wB > 0: B /= wB

    # Stretch is now in units of "bright-star levels", typically ~0.3 - 1.5
    if stretch is None:
        stretch = 0.5

    _vprint(verbose, f"\nLupton parameters: minimum={minimum:.4g}, "
                     f"stretch={stretch:.4g}, Q={Q:.2f}")

    rgb = make_lupton_rgb(R, G, B, minimum=minimum, stretch=stretch, Q=Q)
    if rgb.dtype != np.uint8:
        rgb = (np.clip(rgb, 0, 1) * 255).astype(np.uint8)

    if output_path:
        Image.fromarray(rgb, mode="RGB").save(output_path)
        _vprint(verbose, f"\n✓ Saved colour image to {output_path}")

    if multifilter:
        # In multifilter mode, "filters" is the full list (sorted) and we
        # surface the per-filter weights so callers / PDF subtitles can
        # display them.
        filters_field = active_filters
    else:
        filters_field = list(filter_bands)

    info = {
        "filters": filters_field,
        "rgb_weights": rgb_weights,
        "multifilter": multifilter,
        "wavelengths": {f: wavelengths.get(f) for f in active_filters},
        "n_frames": n_frames,
        "exptimes": exptimes,
        "sky_stats": sky_stats,
        "bg_levels": bg_levels,
        "bright_levels": bright_levels,
        "fwhms_px": fwhms,
        "shifts_px": shifts,
        "pixscale_arcsec": pixscale_arcsec,
        "fov_arcmin": fov_arcmin,
        "center": (center.ra.deg, center.dec.deg),
        "stretch": stretch,
        "Q": Q,
        "minimum": minimum,
        "shape": rgb.shape,
        "method": ("make_lupton_rgb (multi-filter blend)" if multifilter else
                   "make_lupton_rgb (reproject + 2D-bkg + WCS-refine + PSF-match)"),
    }
    return rgb, info


def make_color_from_dataframe(df: pd.DataFrame,
                              center_coord: str,
                              output_dir: str = "/Users/kryanhinds/sedm_phot",
                              cutout_size: int = 1024,
                              seeing_limit: Optional[float] = None,
                              verbose: bool = True) -> Optional[Dict]:
    good = df.copy()
    if seeing_limit is not None and "seeing" in good.columns:
        good = good[(good["seeing"] < seeing_limit) & good["seeing"].notna()].copy()
        _vprint(verbose, f"After seeing filter (<{seeing_limit}\"): "
                          f"{len(good)} / {len(df)} images")
    if len(good) < 3:
        raise ValueError(f"Not enough good images: {len(good)}")
    out = f"{output_dir}/color_composite.png"
    rgb, info = make_color_composite(
        good, center_coord=center_coord, cutout_size=cutout_size,
        output_path=out, verbose=verbose,
    )
    return {"rgb": rgb, "info": info, "df": good, "path": out}


# ---------------------------------------------------------------------------
# Nightly grouping + per-night PDF output
# ---------------------------------------------------------------------------

# Short instrument tag used in PDF filenames. Add to this dict if new
# telescopes are introduced.
INSTRUMENT_TAGS = {
    "NOT/ALFOSC": "NOT",
    "LT/IOO":     "LT",
    "TJO/MEIA2":  "TJO",
    "SEDM":       "SEDM",
    "HCT":        "HCT",
}


def _instrument_tag(name: str) -> str:
    if name in INSTRUMENT_TAGS:
        return INSTRUMENT_TAGS[name]
    # Fallback: take the alphanumeric prefix before the first '/'
    return "".join(c for c in name.split("/")[0] if c.isalnum()) or "UNK"


def find_temporal_groups(df: pd.DataFrame,
                         time_window_days: float = 1.0,
                         instruments: Optional[List[str]] = None,
                         require_multi_instrument: bool = False,
                         min_filters: int = 3) -> List[Tuple[str, pd.DataFrame]]:
    """
    Group frames by a sliding temporal window.

    Two frames are placed in the same group if their ``obs_mjd`` differ by less
    than ``time_window_days``. This is the right tool when you want to combine,
    e.g., a NOT epoch with an LT epoch taken half a night later.

    Parameters
    ----------
    df : DataFrame
        Frame metadata (``filename``, ``filter``, ``instrument``, ``obs_mjd``).
    time_window_days : float
        Half-width of the grouping window in days. ``1.0`` ⇒ ±24 h.
    instruments : list[str], optional
        Restrict to these instruments before grouping.
    require_multi_instrument : bool
        If True, drop groups that contain only one instrument (the original
        ``find_temporal_groups_mjd`` behaviour from the legacy notebook).
    min_filters : int
        Drop groups with fewer than this many distinct filters.

    Returns
    -------
    list[(label, df)]
        ``label`` is a string of the form ``"YYYY-MM-DD"`` derived from the
        median MJD of the group (used for PDF filenames).
    """
    from astropy.time import Time

    d = df.copy()
    if instruments is not None:
        d = d[d["instrument"].isin(instruments)].copy()
    if d.empty:
        return []

    d = d.sort_values("obs_mjd").reset_index(drop=True)
    used = set()
    groups: List[Tuple[str, pd.DataFrame]] = []

    for i, row in d.iterrows():
        if i in used:
            continue
        m = float(row["obs_mjd"])
        mask = (d["obs_mjd"] >= m - time_window_days) & (d["obs_mjd"] <= m + time_window_days)
        g = d[mask].copy()

        if require_multi_instrument and g["instrument"].nunique() < 2:
            continue
        if g["filter"].nunique() < min_filters:
            continue

        used.update(g.index.tolist())
        median_mjd = float(g["obs_mjd"].median())
        label = Time(median_mjd, format="mjd").strftime("%Y-%m-%d")
        groups.append((label, g.reset_index(drop=True)))

    return groups


def find_nightly_groups(df: pd.DataFrame,
                        instruments: Optional[List[str]] = None,
                        min_filters: int = 3) -> List[Tuple[str, pd.DataFrame]]:
    """
    Group rows by night, where a "night" is the integer of ``MJD - 0.5``
    (UTC noon-to-noon). Combines frames from any instrument unless
    ``instruments`` is given.

    Parameters
    ----------
    df : DataFrame
        Must contain ``obs_mjd``, ``filter``, ``instrument``.
    instruments : list[str], optional
        Restrict to these instruments (e.g. ``['NOT/ALFOSC', 'LT/IOO']``).
        Default: all.
    min_filters : int
        Skip nights with fewer than this many distinct filters.

    Returns
    -------
    list[(label, df)]
        ``label`` is a string like ``"2026-03-15"`` (the calendar date of
        the local evening the night belongs to). ``df`` is the frames for
        that night, sorted by MJD.
    """
    from astropy.time import Time

    d = df.copy()
    if instruments is not None:
        d = d[d["instrument"].isin(instruments)].copy()
    if d.empty:
        return []

    d["_night"] = np.floor(d["obs_mjd"].astype(float) - 0.5).astype(int)

    out = []
    for night_mjd, g in d.groupby("_night"):
        if g["filter"].nunique() < min_filters:
            continue
        g = g.sort_values("obs_mjd").reset_index(drop=True)
        # Date of the calendar evening the night belongs to (MJD + 0.5 → UT noon)
        label = Time(int(night_mjd) + 0.5, format="mjd").strftime("%Y-%m-%d")
        out.append((label, g))

    out.sort(key=lambda x: x[0])
    return out


def save_color_pdf(rgb: np.ndarray,
                   info: Dict,
                   output_path: str,
                   target_coord: Optional[str] = None,
                   zoom_arcmin: float = 1.0,
                   circle_arcsec: float = 3.0,
                   title: Optional[str] = None,
                   subtitle: Optional[str] = None,
                   show_scalebar: bool = True,
                   show_compass: bool = True,
                   inset_arcsec: Optional[float] = 20.0,
                   inset_size_frac: float = 0.30,
                   inset_circle_arcsec: Optional[float] = None) -> str:
    """
    Save the RGB image as a PDF cropped to ``zoom_arcmin × zoom_arcmin``
    around the target with a coloured circle marking the transient location.

    Parameters
    ----------
    rgb : (H, W, 3) uint8
        Output of ``make_color_composite``.
    info : dict
        The ``info`` dict returned alongside ``rgb`` (used for pixel scale).
    output_path : str
        Path to the PDF (extension forced to ``.pdf``).
    target_coord : str, optional
        ``"hms dms"`` coordinate string. If omitted, uses the centre of the
        image (which is where ``make_color_composite`` placed its target).
    zoom_arcmin : float
        Side length of the displayed crop, in arcminutes. ``0`` disables
        cropping (full image).
    circle_arcsec : float
        Radius of the marker circle in arcseconds.
    title, subtitle : str, optional
        Title/subtitle text. ``None`` ⇒ no text.
    show_scalebar, show_compass : bool
        Overlay convenience aids.
    inset_arcsec : float, optional
        Side length of a zoom-in inset (in arcsec) drawn inside the top-right
        corner of the main panel, showing the transient + host. ``None`` or
        ``0`` disables the inset. Default ``20`` (i.e. 20"×20").
    inset_size_frac : float
        Width of the inset axes as a fraction of the main panel.
    inset_circle_arcsec : float, optional
        Radius of the marker circle drawn in the inset. Defaults to
        ``circle_arcsec``.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, Rectangle
    from matplotlib.backends.backend_pdf import PdfPages

    pixscale = float(info.get("pixscale_arcsec", 0.3))
    h, w = rgb.shape[:2]

    # Target pixel — make_color_composite always centres the target
    cx_full = (w - 1) / 2.0
    cy_full = (h - 1) / 2.0

    # Crop window
    if zoom_arcmin and zoom_arcmin > 0:
        half_px = zoom_arcmin * 60.0 / 2.0 / pixscale
        x0 = int(max(0, round(cx_full - half_px)))
        x1 = int(min(w, round(cx_full + half_px)))
        y0 = int(max(0, round(cy_full - half_px)))
        y1 = int(min(h, round(cy_full + half_px)))
    else:
        x0, x1, y0, y1 = 0, w, 0, h

    crop = rgb[y0:y1, x0:x1]
    H, W = crop.shape[:2]
    cx, cy = cx_full - x0, cy_full - y0

    if not output_path.lower().endswith(".pdf"):
        output_path = output_path + ".pdf"

    fig, ax = plt.subplots(figsize=(7.0, 7.5))
    ax.imshow(crop, origin="lower", interpolation="nearest")

    # Target circle
    radius_px = max(circle_arcsec / pixscale, 4.0)
    ax.add_patch(Circle((cx, cy), radius=radius_px, fill=False,
                        edgecolor="lime", lw=1.8, alpha=0.95))

    # Scale bar (bottom-left)
    if show_scalebar:
        bar_arcsec = 10.0 if zoom_arcmin >= 0.5 else max(zoom_arcmin * 60.0 / 6.0, 2.0)
        bar_px = bar_arcsec / pixscale
        bx0 = W * 0.06
        by0 = H * 0.06
        ax.plot([bx0, bx0 + bar_px], [by0, by0], color="white", lw=2.2,
                solid_capstyle="butt")
        ax.text(bx0 + bar_px / 2, by0 + H * 0.025,
                f'{int(round(bar_arcsec))}"',
                color="white", ha="center", va="bottom", fontsize=10)

    # N / E compass — placed in the bottom-right (top-right is reserved for the inset,
    # bottom-left for the scale bar, top-left can collide with title text).
    if show_compass:
        L = min(W, H) * 0.06
        ox = W - W * 0.06 - L  # leave room for the E arrow which points LEFT
        oy = H * 0.10
        ax.annotate("", xy=(ox, oy + L), xytext=(ox, oy),
                    arrowprops=dict(arrowstyle="->", color="white", lw=1.4))
        ax.annotate("", xy=(ox - L, oy), xytext=(ox, oy),
                    arrowprops=dict(arrowstyle="->", color="white", lw=1.4))
        ax.text(ox + W * 0.012, oy + L, "N", color="white",
                ha="left", va="center", fontsize=10)
        ax.text(ox - L - W * 0.012, oy, "E", color="white",
                ha="right", va="center", fontsize=10)

    # ------------------------------------------------------------------
    # Zoom-in inset (top-right, inside the main image)
    # ------------------------------------------------------------------
    if inset_arcsec and inset_arcsec > 0:
        # Compute zoom region in main-image pixel coords
        half_inset = inset_arcsec / 2.0 / pixscale
        ix0 = int(max(0, round(cx_full - half_inset)))
        ix1 = int(min(w, round(cx_full + half_inset)))
        iy0 = int(max(0, round(cy_full - half_inset)))
        iy1 = int(min(h, round(cy_full + half_inset)))
        zoom = rgb[iy0:iy1, ix0:ix1]

        # Draw a rectangle on the main panel showing where the inset zooms.
        rect_x0 = ix0 - x0
        rect_y0 = iy0 - y0
        rect_x1 = ix1 - x0
        rect_y1 = iy1 - y0
        if 0 <= rect_x0 < W and 0 <= rect_y0 < H:
            ax.add_patch(Rectangle((rect_x0, rect_y0),
                                   rect_x1 - rect_x0, rect_y1 - rect_y0,
                                   fill=False, edgecolor="white", lw=0.8,
                                   alpha=0.7, ls="--"))

        # Add the inset axes in figure-fraction coordinates anchored to ax
        bbox = ax.get_position()
        pad = 0.018                         # gap from the image edge
        s = float(inset_size_frac)
        ax_in_x = bbox.x1 - s * bbox.width - pad
        ax_in_y = bbox.y1 - s * bbox.height - pad
        ax_in = fig.add_axes([ax_in_x, ax_in_y,
                               s * bbox.width, s * bbox.height])
        ax_in.imshow(zoom, origin="lower", interpolation="nearest")

        # Inset target circle
        inset_radius_arcsec = (inset_circle_arcsec
                               if inset_circle_arcsec is not None
                               else circle_arcsec)
        zh, zw = zoom.shape[:2]
        zcx = cx_full - ix0
        zcy = cy_full - iy0
        zr = max(inset_radius_arcsec / pixscale, 5.0)
        ax_in.add_patch(Circle((zcx, zcy), radius=zr, fill=False,
                               edgecolor="lime", lw=1.8, alpha=0.95))

        # Inset frame (no label)
        for spine in ax_in.spines.values():
            spine.set_edgecolor("white")
            spine.set_linewidth(1.2)
        ax_in.set_xticks([]); ax_in.set_yticks([])

    if title:
        ax.set_title(title, fontsize=12, pad=8)
    if subtitle:
        ax.text(0.5, -0.04, subtitle, transform=ax.transAxes,
                ha="center", va="top", fontsize=9, color="black")

    ax.set_xticks([]); ax.set_yticks([])
    ax.set_aspect("equal")
    with PdfPages(output_path) as pdf:
        pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)
    return output_path


def _render_groups_to_pdfs(groups: List[Tuple[str, pd.DataFrame]],
                           target_coord: str,
                           output_dir: str,
                           target_name: Optional[str],
                           label_suffix: Optional[str],
                           zoom_arcmin: float,
                           cutout_size: int,
                           circle_arcsec: float,
                           inset_arcsec: Optional[float],
                           inset_size_frac: float,
                           verbose: bool,
                           composite_kwargs: Dict) -> List[Dict]:
    """Internal: turn a list of (label, df) groups into one PDF each.

    ``label_suffix``, when set, is appended to the date label in the PDF
    filename. Useful for distinguishing the same date when multiple grouping
    schemes are used in parallel (e.g. ``"window1d"``).
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    results: List[Dict] = []

    for label, g in groups:
        insts = sorted({_instrument_tag(i) for i in g["instrument"].unique()})
        inst_tag = "-".join(insts)

        name_bits = ["color"]
        if target_name:
            name_bits.append(target_name)
        date_part = label if not label_suffix else f"{label}-{label_suffix}"
        name_bits.extend([date_part, inst_tag])
        pdf_path = str(Path(output_dir) / ("_".join(name_bits) + ".pdf"))

        _vprint(verbose, f"\n[{label} | {inst_tag}] {len(g)} frames  "
                          f"filters={sorted(g['filter'].unique())}")
        try:
            rgb, info = make_color_composite(
                g,
                center_coord=target_coord,
                cutout_size=cutout_size,
                output_path=None,
                verbose=verbose,
                **composite_kwargs,
            )
        except Exception as e:
            _vprint(verbose, f"  ✗ composite failed: {e}")
            continue

        n_frames_str = ", ".join(f"{f}:{info['n_frames'][f]}"
                                 for f in info["filters"])
        title = (f"{target_name + '  ·  ' if target_name else ''}"
                 f"{label}  ·  {inst_tag}")
        subtitle = (f"R={info['filters'][0]}  G={info['filters'][1]}  "
                    f"B={info['filters'][2]}   |   "
                    f"frames stacked: {n_frames_str}   |   "
                    f"crop: {zoom_arcmin:.2f}′")

        save_color_pdf(
            rgb, info, pdf_path,
            target_coord=target_coord,
            zoom_arcmin=zoom_arcmin,
            circle_arcsec=circle_arcsec,
            inset_arcsec=inset_arcsec,
            inset_size_frac=inset_size_frac,
            title=title,
            subtitle=subtitle,
        )
        _vprint(verbose, f"  ✓ saved {pdf_path}")
        results.append({
            "label": label,
            "instruments": insts,
            "n_frames": int(g.shape[0]),
            "filters": info["filters"],
            "path": pdf_path,
        })
    return results


def make_nightly_pdfs(df: pd.DataFrame,
                      target_coord: str,
                      output_dir: str,
                      instruments: Optional[List[str]] = None,
                      target_name: Optional[str] = None,
                      zoom_arcmin: float = 1.0,
                      cutout_size: int = 1024,
                      circle_arcsec: float = 3.0,
                      inset_arcsec: Optional[float] = 20.0,
                      inset_size_frac: float = 0.30,
                      min_filters: int = 3,
                      verbose: bool = True,
                      **composite_kwargs) -> List[Dict]:
    """
    Build one nightly PDF per night that has ``min_filters`` filters.

    Filename pattern
    ----------------
    ``<output_dir>/color_<TARGET>_<DATE>_<INSTRUMENTS>.pdf`` -- e.g.
    ``color_ZTF26aakjzdt_2026-03-15_LT-NOT.pdf``. ``TARGET`` is omitted if
    ``target_name`` is None.

    Parameters
    ----------
    df : DataFrame
        Frame metadata (``filename``, ``filter``, ``instrument``, ``obs_mjd``,
        ...) -- typically produced by the loader cell of the notebook.
    target_coord : str
        ``"hms dms"`` coord string passed to :func:`make_color_composite`.
    output_dir : str
        Directory the PDFs are saved into.
    instruments : list[str], optional
        Restrict to these instruments. Default: all.
    target_name : str, optional
        If set, included in the filename.
    zoom_arcmin : float
        Display crop size in arcmin.
    cutout_size : int
        Internal pipeline cutout size in pixels (full reduction grid).
    circle_arcsec : float
        Marker radius in arcseconds.
    min_filters : int
        Skip nights with fewer than this many filters.
    composite_kwargs
        Forwarded to :func:`make_color_composite` (``Q``, ``stretch``,
        ``channel_weights`` …).

    Returns
    -------
    list[dict]
        One entry per produced PDF: ``{date, instruments, n_frames, path}``.
    """
    nights = find_nightly_groups(df, instruments=instruments,
                                 min_filters=min_filters)
    _vprint(verbose, f"Nightly groups (≥{min_filters} filters): {len(nights)}")
    return _render_groups_to_pdfs(
        nights,
        target_coord=target_coord,
        output_dir=output_dir,
        target_name=target_name,
        label_suffix=None,
        zoom_arcmin=zoom_arcmin,
        cutout_size=cutout_size,
        circle_arcsec=circle_arcsec,
        inset_arcsec=inset_arcsec,
        inset_size_frac=inset_size_frac,
        verbose=verbose,
        composite_kwargs=composite_kwargs,
    )


def make_temporal_group_pdfs(df: pd.DataFrame,
                             target_coord: str,
                             output_dir: str,
                             time_window_days: float = 1.0,
                             instruments: Optional[List[str]] = None,
                             require_multi_instrument: bool = False,
                             target_name: Optional[str] = None,
                             zoom_arcmin: float = 1.0,
                             cutout_size: int = 1024,
                             circle_arcsec: float = 3.0,
                             inset_arcsec: Optional[float] = 20.0,
                             inset_size_frac: float = 0.30,
                             min_filters: int = 3,
                             label_suffix: Optional[str] = None,
                             verbose: bool = True,
                             **composite_kwargs) -> List[Dict]:
    """
    Build one PDF per temporal group (sliding ±``time_window_days`` window).

    This is the original ``find_temporal_groups_mjd`` behaviour from the legacy
    notebook, wrapped as a one-call function. Use this when you want to combine
    frames taken within a few hours / a day of each other across multiple
    instruments (e.g., a NOT epoch and an LT epoch from the same observing
    campaign that didn't quite overlap on the same calendar night).

    Parameters
    ----------
    df : DataFrame
        Frame metadata.
    target_coord : str
        ``"hms dms"`` coordinate string.
    output_dir : str
        Directory the PDFs are saved to.
    time_window_days : float
        Half-width of the grouping window. ``1.0`` ⇒ ±24 h.
    instruments : list[str], optional
        Restrict to these instruments before grouping.
    require_multi_instrument : bool
        If True, drop groups that contain only one instrument (matches the
        legacy behaviour where you only wanted multi-telescope epochs).
    target_name, zoom_arcmin, cutout_size, circle_arcsec, min_filters,
    composite_kwargs
        Same as :func:`make_nightly_pdfs`.
    label_suffix : str, optional
        Appended to the filename's date part (default ``"win{N}d"`` if ``N``
        is non-default). Useful when you produce both nightly and temporal
        PDFs into the same folder and want to keep them distinct.

    Returns
    -------
    list[dict]
    """
    groups = find_temporal_groups(
        df,
        time_window_days=time_window_days,
        instruments=instruments,
        require_multi_instrument=require_multi_instrument,
        min_filters=min_filters,
    )
    _vprint(verbose, f"Temporal groups (±{time_window_days} d, "
                     f"≥{min_filters} filters"
                     f"{', multi-instrument only' if require_multi_instrument else ''}"
                     f"): {len(groups)}")

    if label_suffix is None and (
        time_window_days != 1.0 or require_multi_instrument):
        bits = [f"win{time_window_days:g}d"]
        if require_multi_instrument:
            bits.append("multi")
        label_suffix = "-".join(bits)

    return _render_groups_to_pdfs(
        groups,
        target_coord=target_coord,
        output_dir=output_dir,
        target_name=target_name,
        label_suffix=label_suffix,
        zoom_arcmin=zoom_arcmin,
        cutout_size=cutout_size,
        circle_arcsec=circle_arcsec,
        inset_arcsec=inset_arcsec,
        inset_size_frac=inset_size_frac,
        verbose=verbose,
        composite_kwargs=composite_kwargs,
    )


# ---------------------------------------------------------------------------
# Multi-panel finder figure (colour + per-filter sci/ref/sub cutouts)
# ---------------------------------------------------------------------------

def _filter_instrument_breakdown(df: pd.DataFrame) -> str:
    """Return a "INST1: filters, INST2: filters" string, e.g.
    ``'NOT/ALFOSC: griz, LT/IOO: gri'``.
    """
    if df is None or len(df) == 0 or "instrument" not in df.columns:
        return ""
    parts = []
    # Sort instruments by typical wavelength of the bluest filter (just for
    # visual stability — astronomy convention is bluest first).
    insts = sorted(df["instrument"].dropna().unique())
    for inst in insts:
        grp = df[df["instrument"] == inst]
        flts = sorted({_normalise_filter(f) for f in grp["filter"].dropna()
                       if _normalise_filter(f)},
                      key=lambda f: DEFAULT_WAVELENGTHS.get(f, 0))
        if flts:
            parts.append(f"{inst}: {''.join(flts)}")
    return ", ".join(parts)


def _date_str_from_df(df: pd.DataFrame) -> str:
    """Median-MJD → 'YYYY-MM-DD'. Empty string if no time info."""
    if df is None or "obs_mjd" not in df.columns:
        return ""
    mjd = df["obs_mjd"].dropna()
    mjd = mjd[mjd > 0]
    if not len(mjd):
        return ""
    from astropy.time import Time
    return Time(float(mjd.median()), format="mjd").strftime("%Y-%m-%d")


def _draw_color_into_axes(fig, ax, rgb, info,
                          target_coord: str,
                          zoom_arcmin: float,
                          inset_arcsec: Optional[float],
                          inset_size_frac: float,
                          circle_arcsec: float,
                          label: Optional[str] = None) -> None:
    """Draw a colour composite (with target circle + optional zoom inset) into
    the given axes, square aspect."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, Rectangle

    pixscale = float(info.get("pixscale_arcsec", 0.3))
    h, w = rgb.shape[:2]
    cx_full = (w - 1) / 2.0
    cy_full = (h - 1) / 2.0

    # Crop
    if zoom_arcmin and zoom_arcmin > 0:
        half_px = zoom_arcmin * 60.0 / 2.0 / pixscale
        x0 = int(max(0, round(cx_full - half_px)))
        x1 = int(min(w, round(cx_full + half_px)))
        y0 = int(max(0, round(cy_full - half_px)))
        y1 = int(min(h, round(cy_full + half_px)))
    else:
        x0, x1, y0, y1 = 0, w, 0, h
    crop = rgb[y0:y1, x0:x1]
    H, W = crop.shape[:2]

    ax.imshow(crop, origin="lower", interpolation="nearest")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_linewidth(2.0)

    cx, cy = cx_full - x0, cy_full - y0
    rpx = max(circle_arcsec / pixscale, 4.0)
    ax.add_patch(Circle((cx, cy), radius=rpx, fill=False,
                        edgecolor="lime", lw=1.8, alpha=0.95))

    if inset_arcsec and inset_arcsec > 0:
        half_inset = inset_arcsec / 2.0 / pixscale
        ix0 = int(max(0, round(cx_full - half_inset)))
        ix1 = int(min(w, round(cx_full + half_inset)))
        iy0 = int(max(0, round(cy_full - half_inset)))
        iy1 = int(min(h, round(cy_full + half_inset)))
        zoom = rgb[iy0:iy1, ix0:ix1]

        # Dashed rectangle on the main panel showing where the inset zooms.
        rect_x0 = ix0 - x0
        rect_y0 = iy0 - y0
        rect_x1 = ix1 - x0
        rect_y1 = iy1 - y0
        if 0 <= rect_x0 < W and 0 <= rect_y0 < H:
            ax.add_patch(Rectangle((rect_x0, rect_y0),
                                   rect_x1 - rect_x0, rect_y1 - rect_y0,
                                   fill=False, edgecolor="white", lw=0.8,
                                   alpha=0.7, ls="--"))

        # Use Axes.inset_axes() in axes-fraction coords. This places the
        # inset in the top-right corner of the parent axes regardless of
        # any subsequent layout adjustments — the previous fig.add_axes()
        # approach was sensitive to draw-order and figure-layout state.
        s = float(inset_size_frac)
        ax_in = ax.inset_axes([1.0 - s, 1.0 - s, s, s])
        ax_in.imshow(zoom, origin="lower", interpolation="nearest")
        ax_in.set_aspect("equal", adjustable="box")
        ax_in.set_xticks([]); ax_in.set_yticks([])
        for spine in ax_in.spines.values():
            spine.set_edgecolor("white")
            spine.set_linewidth(1.2)

        zcx = cx_full - ix0
        zcy = cy_full - iy0
        zr = max(circle_arcsec / pixscale, 5.0)
        ax_in.add_patch(Circle((zcx, zcy), radius=zr, fill=False,
                               edgecolor="lime", lw=1.8, alpha=0.95))

    if label:
        ax.text(0.5, 0.025, label, transform=ax.transAxes,
                ha="center", va="bottom", fontsize=14, color="white",
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.4",
                          facecolor="black", edgecolor="none", alpha=0.7))


def make_finder_figure(science_df: pd.DataFrame,
                        ref_df: pd.DataFrame,
                        cutouts: Dict[str, Dict[str, str]],
                        target_coord: str,
                        output_path: str,
                        filters: Optional[List[str]] = None,
                        rows: Optional[List[Tuple[str, str]]] = None,
                        target_name: Optional[str] = None,
                        science_label: Optional[str] = None,
                        ref_label: Optional[str] = None,
                        color_zoom_arcmin: float = 1.0,
                        color_inset_arcsec: Optional[float] = 20.0,
                        color_inset_size_frac: float = 0.30,
                        color_cutout_size_sci: Optional[int] = None,
                        color_cutout_size_ref: Optional[int] = None,
                        sci_filter_bands=None,
                        ref_filter_bands="all",
                        circle_arcsec: float = 3.0,
                        cutout_circle_arcsec: float = 2.0,
                        cutout_scale_arcsec: float = 5.0,
                        cutout_radius_pix: int = 50,
                        top_panel_size: float = 5.5,
                        bottom_panel_size: Optional[float] = None,
                        cmap: str = "viridis",
                        figsize: Optional[Tuple[float, float]] = None,
                        verbose: bool = True,
                        **composite_kwargs) -> str:
    """
    Make a finder-style PDF that shows, in one figure:

      • TOP — two large square colour-composite panels (reference left,
        science right), each with a target circle and a zoom inset.
      • BOTTOM — three rows × N columns of square greyscale cutouts
        (Science / Reference / Difference) for every filter listed in
        ``cutouts``, with a circle marking the transient location.

    Parameters
    ----------
    science_df, ref_df : DataFrame
        Frame metadata for the top-row colour composites (same shape as the
        DataFrames you'd pass to ``make_color_composite``).
    cutouts : dict
        ``{'g': {'sci': path, 'ref': path, 'sub': path}, 'r': {...}, ...}``.
        Filenames pointing to per-filter cutout FITS at the transient
        location.
    filters : list[str], optional
        Which filters to show at the bottom and in what column order, e.g.
        ``['r', 'i', 'z']`` or ``['g', 'r', 'i', 'z']``. Filters not in
        this list are dropped from the figure even if present in
        ``cutouts``. Default ``None`` ⇒ use every filter in ``cutouts`` in
        its dict order.
    rows : list[(label, key)], optional
        Which rows to render at the bottom and in what order. Each tuple
        is ``(display_label, dict_key)`` where ``dict_key`` matches a key
        inside ``cutouts[filter]``. Default
        ``[('Science','sci'), ('Reference','ref'), ('Difference','sub')]``.
        Pass e.g. ``[('Science','sci'), ('Difference','sub')]`` to drop
        the reference row.
    target_coord : str
        ``"hms dms"`` coordinate string. Used to centre the colour panels and
        to place the target circle in each cutout (uses each cutout's own
        WCS, falling back to the image centre if the WCS isn't usable).
    output_path : str
        PDF path.
    science_label, ref_label : str, optional
        Tag drawn in the top-left of each colour panel (e.g. ``"T₀+2.1 d"``,
        ``"PS1 + LS reference"``).
    color_zoom_arcmin : float
        Side length of the displayed colour panels **on the sky**, in
        arcminutes. e.g. ``1.0`` ⇒ 60″×60″ field shown in each panel.
    color_inset_arcsec : float, optional
        Side length of the **zoom inset** (top-right of each colour panel)
        on the sky, in arcseconds. Set ``None`` or ``0`` to disable.
    color_inset_size_frac : float
        Width of the inset axes as a fraction of the parent colour panel.
    color_cutout_size_sci, color_cutout_size_ref : int, optional
        Internal pipeline grid sizes in **pixels** for the colour composites.
        These do NOT control the displayed image size — that's
        ``color_zoom_arcmin``. They control how much sky the pipeline reduces
        before cropping. Default ``None`` ⇒ auto-sized so the pipeline grid
        is roughly 3× ``color_zoom_arcmin`` (always large enough to cover
        the visible crop).
    sci_filter_bands, ref_filter_bands
        Forwarded to ``make_color_composite``. Default for refs is ``"all"``
        (multi-filter wavelength blend); science uses the auto-pick of the
        three reddest filters.
    cutout_circle_arcsec, cutout_scale_arcsec : float
        Marker / scale-bar sizes for the bottom cutouts.
    cutout_radius_pix : int
        Half-side of the bottom-panel cutouts in pixels (uses each FITS file's
        own WCS to centre on the transient). Default ``50`` ⇒ each panel shows
        a 100×100 px square around the target.
    top_panel_size : float
        Side length (inches) of each of the two square colour panels.
        Default ``5.5``.
    bottom_panel_size : float, optional
        Side length (inches) of each greyscale cutout cell. Default
        ``2 * top_panel_size / nf`` so the colour panels and cutout row span
        the same total width and every cell is square. Set explicitly (e.g.
        ``2.5``) to make the cutouts smaller; the cutout row will then sit
        centred under the wider colour-panel row.
    cmap : str
        Matplotlib colour map for the greyscale cutouts (default ``"viridis"``
        to match your existing finder style).
    figsize : (W, H), optional
        Defaults to ``(3*ncols, 3*nrows + ...)`` so all cells are square.
    composite_kwargs
        Forwarded to ``make_color_composite`` (``Q``, ``stretch`` …).

    Returns
    -------
    str
        The path the PDF was written to.
    """
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from matplotlib.patches import Circle
    from astropy.coordinates import SkyCoord
    from astropy.visualization import ZScaleInterval

    target = SkyCoord(target_coord, unit=("hourangle", "deg"))

    # Filter composite_kwargs to those make_color_composite actually accepts,
    # so users can sprinkle finder-only kwargs into the call without crashing.
    import inspect
    valid = set(inspect.signature(make_color_composite).parameters.keys())
    dropped = [k for k in composite_kwargs.keys() if k not in valid]
    if dropped:
        _vprint(verbose, f"⚠ ignoring unknown kwargs (not for make_color_composite): "
                          f"{dropped}")
    composite_kwargs = {k: v for k, v in composite_kwargs.items() if k in valid}

    # Resolve filter order
    if filters is None:
        filters = list(cutouts.keys())
    else:
        filters = list(filters)
        missing = [f for f in filters if f not in cutouts]
        if missing:
            _vprint(verbose, f"⚠ filters {missing} not in cutouts dict — they "
                              f"will render as empty panels")
    if not filters:
        raise ValueError("`filters` resolved to an empty list "
                         "(empty cutouts and no filters arg)")
    nf = len(filters)
    half = max(nf // 2, 1)

    # Resolve rows
    if rows is None:
        rows = [("Science", "sci"), ("Reference", "ref"), ("Difference", "sub")]
    rows = list(rows)
    nrows_bottom = len(rows)

    # ------------------------------------------------------------------
    # 1. Build colour composites for the top row
    #
    # color_cutout_size_* is the internal pipeline grid in pixels (NOT the
    # displayed crop -- that's color_zoom_arcmin). If unspecified we pick a
    # size that gives ~3× headroom around the zoom_arcmin crop, so the
    # composite always covers the visible area with margin.
    # ------------------------------------------------------------------
    def _auto_pipe_size(default_size: int, ps_guess_arcsec: float = 0.25) -> int:
        if color_zoom_arcmin and color_zoom_arcmin > 0:
            need = int(np.ceil(color_zoom_arcmin * 60.0 / ps_guess_arcsec * 3.0))
            return max(default_size, need)
        return default_size

    _ref_size = (color_cutout_size_ref if color_cutout_size_ref is not None
                 else _auto_pipe_size(600))
    _sci_size = (color_cutout_size_sci if color_cutout_size_sci is not None
                 else _auto_pipe_size(1024, ps_guess_arcsec=0.21))

    # Reference cutouts (PS1, Legacy Survey, ...) are already coadds and
    # essentially CR-free; running LACosmic on them is slow and can clip
    # real star cores. Default the reference call to lacosmic=False unless
    # the user has explicitly opted in via composite_kwargs.
    ref_kwargs = dict(composite_kwargs)
    ref_kwargs.setdefault("lacosmic", False)

    _vprint(verbose, f"Building reference colour composite "
                     f"(internal grid {_ref_size}×{_ref_size}px)...")
    rgb_ref, info_ref = make_color_composite(
        ref_df, filter_bands=ref_filter_bands,
        center_coord=target_coord, cutout_size=_ref_size,
        verbose=False, **ref_kwargs,
    )
    _vprint(verbose, f"Building science colour composite "
                     f"(internal grid {_sci_size}×{_sci_size}px)...")
    rgb_sci, info_sci = make_color_composite(
        science_df, filter_bands=sci_filter_bands,
        center_coord=target_coord, cutout_size=_sci_size,
        verbose=False, **composite_kwargs,
    )

    # Sanity check: warn if the user-requested zoom exceeds the pipeline grid
    for tag, info in (("reference", info_ref), ("science", info_sci)):
        fov_pipe = info.get("fov_arcmin", 0)
        if color_zoom_arcmin > fov_pipe:
            _vprint(verbose, f"⚠ {tag} colour panel: requested zoom "
                              f"({color_zoom_arcmin:.2f}′) exceeds the pipeline "
                              f"grid FOV ({fov_pipe:.2f}′) — the panel will be "
                              f"clamped to the available data. Increase "
                              f"color_cutout_size_{tag[:3]} to fix.")

    # ------------------------------------------------------------------
    # 2. Layout — two INDEPENDENT GridSpecs stacked vertically:
    #      top: 1 row × 2 cols (the two square colour panels)
    #      bot: nrows_bottom × nf (the smaller square cutout cells)
    #
    #    Decoupling means the cutout aspect / data shape can never throw off
    #    the colour-panel sizing.
    # ------------------------------------------------------------------
    if bottom_panel_size is None:
        # Default: bottom row shares the same total width as the top row,
        # giving square cutout cells with bottom_size = 2*top_size / nf.
        bottom_panel_size = 2.0 * top_panel_size / nf

    top_section_w = 2.0 * top_panel_size
    bot_section_w = nf * bottom_panel_size
    fig_width = max(top_section_w, bot_section_w)
    top_section_h = top_panel_size
    bot_section_h = nrows_bottom * bottom_panel_size

    if figsize is None:
        figsize = (fig_width, top_section_h + bot_section_h)

    fig = plt.figure(figsize=figsize)

    # The two sections live in figure-fraction coords. We use add_gridspec
    # with explicit left/right so the bot can be narrower than top (or vice
    # versa) and still appear centred.
    actual_w, actual_h = fig.get_size_inches()
    top_left = (actual_w - top_section_w) / 2.0 / actual_w
    top_right = top_left + top_section_w / actual_w
    bot_left = (actual_w - bot_section_w) / 2.0 / actual_w
    bot_right = bot_left + bot_section_w / actual_w
    bot_top = bot_section_h / actual_h
    bot_bottom = 0.0
    top_bottom = bot_top
    top_top = 1.0

    top_gs = fig.add_gridspec(
        1, 2, left=top_left, right=top_right,
        bottom=top_bottom, top=top_top, wspace=0.0,
    )
    bot_gs = fig.add_gridspec(
        nrows_bottom, nf, left=bot_left, right=bot_right,
        bottom=bot_bottom, top=bot_top, hspace=0.0, wspace=0.0,
    )

    ax_ref_color = fig.add_subplot(top_gs[0, 0])
    ax_sci_color = fig.add_subplot(top_gs[0, 1])

    # Auto-build the labels with per-instrument filter breakdown if the user
    # didn't supply one. Format examples:
    #   "Reference  ·  PS1: grizy, LegacySurvey: grz"
    #   "Science  ·  2026-04-23  ·  NOT/ALFOSC: griz, LT/IOO: gri"
    if ref_label is None:
        breakdown = _filter_instrument_breakdown(ref_df)
        ref_label = f"Reference  ·  {breakdown}" if breakdown else "Reference"
    if science_label is None:
        breakdown = _filter_instrument_breakdown(science_df)
        date_part = _date_str_from_df(science_df)
        # Line 1: "Science  ·  2026-04-23"
        line1_parts = ["Science"]
        if date_part:
            line1_parts.append(date_part)
        line1 = "  ·  ".join(line1_parts)
        # Line 2: per-instrument filter breakdown
        science_label = f"{line1}\n{breakdown}" if breakdown else line1

    _draw_color_into_axes(
        fig, ax_ref_color, rgb_ref, info_ref,
        target_coord=target_coord,
        zoom_arcmin=color_zoom_arcmin,
        inset_arcsec=color_inset_arcsec,
        inset_size_frac=color_inset_size_frac,
        circle_arcsec=circle_arcsec,
        label=ref_label,
    )
    _draw_color_into_axes(
        fig, ax_sci_color, rgb_sci, info_sci,
        target_coord=target_coord,
        zoom_arcmin=color_zoom_arcmin,
        inset_arcsec=color_inset_arcsec,
        inset_size_frac=color_inset_size_frac,
        circle_arcsec=circle_arcsec,
        label=science_label,
    )

    # ------------------------------------------------------------------
    # 3. Bottom rows × nf cols — greyscale cutouts driven by `rows`
    # ------------------------------------------------------------------
    for col, filt in enumerate(filters):
        per_filter = cutouts.get(filt, {})
        for row, (label, key) in enumerate(rows):
            ax = fig.add_subplot(bot_gs[row, col])
            ax.set_aspect("equal", adjustable="box")
            ax.set_xticks([]); ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_linewidth(2.0)

            path = per_filter.get(key)
            data = wcs = None
            if path:
                try:
                    data, wcs, _hdr = _open_image(path)
                except Exception as e:
                    _vprint(verbose, f"  ✗ {filt}/{key} ({path}): {e}")
                    data, wcs = None, None

            if data is not None:
                # Crop to ±cutout_radius_pix around the target so every panel
                # is a small square zoomed onto the transient + host. Without
                # this, the FITS files (often the entire field) drown the
                # transient in stars and break the figure layout.
                cx_target = (data.shape[1] - 1) / 2.0
                cy_target = (data.shape[0] - 1) / 2.0
                if wcs is not None:
                    try:
                        x_t, y_t = wcs.world_to_pixel(target)
                        cx_target, cy_target = float(x_t), float(y_t)
                    except Exception:
                        pass

                ix0 = int(round(cx_target - cutout_radius_pix))
                ix1 = int(round(cx_target + cutout_radius_pix))
                iy0 = int(round(cy_target - cutout_radius_pix))
                iy1 = int(round(cy_target + cutout_radius_pix))
                ix0 = max(0, ix0); ix1 = min(data.shape[1], ix1)
                iy0 = max(0, iy0); iy1 = min(data.shape[0], iy1)
                if ix1 - ix0 >= 4 and iy1 - iy0 >= 4:
                    data = data[iy0:iy1, ix0:ix1]
                    cx_px = cx_target - ix0
                    cy_px = cy_target - iy0
                else:
                    cx_px, cy_px = cx_target, cy_target

                finite = data[np.isfinite(data)]
                is_diff = key in ("sub", "diff", "difference")
                if is_diff:
                    vlim = float(np.nanpercentile(np.abs(finite), 98)) if finite.size else 1.0
                    ax.imshow(data, origin="lower", cmap=cmap,
                              vmin=-vlim, vmax=vlim)
                else:
                    if finite.size:
                        vmin, vmax = ZScaleInterval().get_limits(finite)
                    else:
                        vmin, vmax = 0.0, 1.0
                    ax.imshow(data, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)

                # Target circle
                rpx = 6.0
                if wcs is not None:
                    try:
                        ps = _pixel_scale_arcsec(wcs)
                        rpx = max(cutout_circle_arcsec / ps, 4.0)
                    except Exception:
                        pass
                ax.add_patch(Circle((cx_px, cy_px), radius=rpx, fill=False,
                                    edgecolor="magenta", lw=1.5, alpha=0.95))

                # Scale bar on top row only
                if row == 0 and wcs is not None:
                    try:
                        ps = _pixel_scale_arcsec(wcs)
                        bar_px = cutout_scale_arcsec / ps
                        ny, nx = data.shape
                        bx0 = nx * 0.06
                        by0 = ny * 0.08
                        ax.plot([bx0, bx0 + bar_px], [by0, by0], "w-", lw=2.5,
                                solid_capstyle="butt")
                        ax.text(bx0 + bar_px / 2, by0 + ny * 0.025,
                                f'{int(round(cutout_scale_arcsec))}"',
                                color="white", ha="center", va="bottom",
                                fontsize=10)
                    except Exception:
                        pass

            # Labels
            if col == 0:
                ax.text(0.04, 0.96, label, transform=ax.transAxes,
                        va="top", ha="left", fontsize=11, color="white",
                        bbox=dict(boxstyle="round,pad=0.25",
                                  facecolor="black", edgecolor="none", alpha=0.6))
            if row == 0:
                ax.text(0.96, 0.96, f'{filt}-band', transform=ax.transAxes,
                        va="top", ha="right", fontsize=10, color="navy",
                        fontweight="bold",
                        bbox=dict(boxstyle="round,pad=0.25",
                                  facecolor="white", edgecolor="none", alpha=0.85))

    if not output_path.lower().endswith(".pdf"):
        output_path = output_path + ".pdf"
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    _vprint(verbose, f"\n✓ Saved {output_path}")
    return output_path
