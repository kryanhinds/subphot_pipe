"""
build_psf.py — Robust empirical ePSF construction using photutils EPSFBuilder.

Replaces the simple median-stacking approach in psf_measure.py with the
Anderson & King (2000, 2006) iterative algorithm.  Key improvements:

  - sep-based 2-D background subtraction (handles gradients and bright haloes)
  - Isolation filtering: stars with neighbours within isolation_factor × FWHM
    are excluded before the ePSF is built
  - Iterative ePSF construction (EPSFBuilder): centroid offsets and PSF shape
    are refined together across iterations, eliminating the dipole residuals
    that arise from simple median-stacking of mis-centred cutouts
  - Residual chi² quality metric directly comparable to PSFEx CHI2
  - Native pixel-scale kernel output ready for scipy.signal.fftconvolve

Pipeline drop-in usage
----------------------
from build_psf import build_psf_from_fits

result = build_psf_from_fits('aligned_science.fits', fwhm_guess=4.5)

if result['chi2'] < 10:                  # quality gate (PSFEx uses ~3)
    kernel = result['kernel']            # 2-D array, normalised, native px
else:
    # fall back to Gaussian
    ...

Returns
-------
dict with keys:
  'kernel'     : 2-D ndarray  — flux-normalised PSF at native pixel scale
  'epsf'       : EPSFModel    — photutils model (oversampled, for inspection)
  'fwhm'       : float        — circularised FWHM = sqrt(fwhm_x * fwhm_y) [px]
  'fwhm_x'     : float        — FWHM along x [px]
  'fwhm_y'     : float        — FWHM along y [px]
  'theta'      : float        — position angle of PSF major axis [rad]
  'elongation' : float        — fwhm_major / fwhm_minor  (1 = round)
  'chi2'       : float        — residual chi² of ePSF fit (1 = perfect)
  'n_stars'    : int          — number of stars used
  'stars_x'    : ndarray      — x pixel positions of used stars
  'stars_y'    : ndarray      — y pixel positions of used stars
  'oversampling': int         — oversampling factor used in EPSFBuilder
"""

import warnings
import numpy as np
from astropy.io import fits
from astropy.nddata import NDData
from astropy.stats import sigma_clipped_stats, SigmaClip
from astropy.table import Table
from scipy.optimize import curve_fit
from scipy.ndimage import zoom

try:
    import sep
    _HAS_SEP = True
except ImportError:
    _HAS_SEP = False

from photutils.detection import IRAFStarFinder
from photutils.psf import EPSFBuilder, extract_stars
from photutils.background import Background2D, SExtractorBackground


# ─────────────────────────────────────────────────────────────────────────────
# Background estimation
# ─────────────────────────────────────────────────────────────────────────────

def _estimate_background(data, box_size=64, filter_size=3):
    """
    Estimate a smooth 2-D background using SExtractorBackground via
    photutils.  Falls back to sep if photutils raises, then to sigma-clipped
    global median if both fail.

    Returns (bkg_2d, rms_2d) — both same shape as data.
    """
    ny, nx = data.shape
    box = max(32, min(box_size, ny // 6, nx // 6))
    # Ensure box divides reasonably
    box = int(box)
    filt = max(1, min(filter_size, box // 8))

    # ── photutils Background2D ────────────────────────────────────────────────
    try:
        sigma_clip = SigmaClip(sigma=3.0, maxiters=5)
        bkg = Background2D(
            data, box_size=(box, box),
            filter_size=(filt, filt),
            sigma_clip=sigma_clip,
            bkg_estimator=SExtractorBackground(),
        )
        return bkg.background, bkg.background_rms
    except Exception:
        pass

    # ── sep fallback ──────────────────────────────────────────────────────────
    if _HAS_SEP:
        try:
            d = np.ascontiguousarray(data, dtype=np.float64)
            bkg_sep = sep.Background(d, bw=box, bh=box, fw=filt, fh=filt)
            bkg_arr = bkg_sep.back()
            rms_arr = bkg_sep.rms()
            return bkg_arr, rms_arr
        except Exception:
            pass

    # ── global sigma-clipped median fallback ──────────────────────────────────
    _, med, rms = sigma_clipped_stats(data, sigma=3.0, maxiters=5)
    return np.full_like(data, med, dtype=float), np.full_like(data, rms, dtype=float)


# ─────────────────────────────────────────────────────────────────────────────
# Source detection & isolation filtering
# ─────────────────────────────────────────────────────────────────────────────

def _detect_stars(data_sub, rms, fwhm_guess, threshold_sigma=5.0,
                  saturation=None, edge_margin=None):
    """
    Detect point sources with IRAFStarFinder and apply basic quality cuts.

    Returns an astropy Table with columns xcentroid, ycentroid, peak, flux.
    """
    ny, nx = data_sub.shape
    margin = int(edge_margin if edge_margin is not None else max(20, fwhm_guess * 4))
    threshold = threshold_sigma * float(np.median(rms))

    finder = IRAFStarFinder(
        threshold=threshold,
        fwhm=fwhm_guess,
        roundhi=0.5,          # reject very elongated (galaxies, cosmic rays)
        sharplo=0.3,
        sharphi=0.9,
        peakmax=saturation,   # IRAFStarFinder ignores None automatically
    )

    # Trim edges before detection so the finder doesn't find edge artefacts
    inner = data_sub[margin: ny - margin, margin: nx - margin]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        sources = finder(inner)

    if sources is None or len(sources) == 0:
        return Table(names=['xcentroid', 'ycentroid', 'peak', 'flux'],
                     dtype=[float, float, float, float])

    # Offset coordinates back to full-frame
    sources['xcentroid'] += margin
    sources['ycentroid'] += margin

    # Remove saturated stars
    if saturation is not None:
        sources = sources[sources['peak'] < saturation]

    return sources


def _filter_isolated(sources, fwhm, isolation_factor=2.5):
    """
    Keep only stars that have no neighbour within isolation_factor * fwhm pixels.
    This is critical for ePSF construction — contaminated cutouts bias the PSF.
    """
    if len(sources) < 2:
        return sources

    x = np.array(sources['xcentroid'])
    y = np.array(sources['ycentroid'])
    min_sep = isolation_factor * fwhm

    keep = np.ones(len(x), dtype=bool)
    for i in range(len(x)):
        for j in range(len(x)):
            if i == j:
                continue
            if np.hypot(x[i] - x[j], y[i] - y[j]) < min_sep:
                keep[i] = False
                break

    return sources[keep]


# ─────────────────────────────────────────────────────────────────────────────
# 2-D Gaussian fit to the final kernel (for FWHM measurement)
# ─────────────────────────────────────────────────────────────────────────────

def _gaussian_2d(xy, amplitude, x0, y0, sigma_x, sigma_y, theta, offset):
    x, y = xy
    dx, dy = x - x0, y - y0
    ct, st = np.cos(theta), np.sin(theta)
    a = (ct / sigma_x) ** 2 + (st / sigma_y) ** 2
    b = np.sin(2 * theta) * (1 / sigma_x ** 2 - 1 / sigma_y ** 2)
    c = (st / sigma_x) ** 2 + (ct / sigma_y) ** 2
    return (offset + amplitude * np.exp(-0.5 * (a * dx**2 + b * dx * dy + c * dy**2))).ravel()


def _fit_psf_fwhm(kernel):
    """
    Fit a 2-D tilted Gaussian to the PSF kernel.
    Returns (fwhm_x, fwhm_y, theta_rad).
    """
    ny, nx = kernel.shape
    y_g, x_g = np.mgrid[:ny, :nx]
    total = kernel.sum()
    if total <= 0:
        return float(nx / 4), float(ny / 4), 0.0

    cx = float((x_g * kernel).sum() / total)
    cy = float((y_g * kernel).sum() / total)
    dx = x_g - cx
    dy = y_g - cy
    sx = float(np.sqrt(max(0.5, (dx**2 * kernel).sum() / total)))
    sy = float(np.sqrt(max(0.5, (dy**2 * kernel).sum() / total)))
    amp = float(kernel.max())
    p0 = [amp, cx, cy, sx, sy, 0.0, 0.0]
    lim = min(nx, ny)
    bounds = ([0, -1, -1, 0.3, 0.3, -np.pi/2, -np.inf],
              [np.inf, nx+1, ny+1, lim, lim, np.pi/2, np.inf])
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            popt, _ = curve_fit(_gaussian_2d, (x_g, y_g), kernel.ravel(),
                                p0=p0, bounds=bounds, maxfev=5000,
                                ftol=1e-7, xtol=1e-7)
        _, _, _, sx_f, sy_f, theta_f, _ = popt
        fwhm_x = 2.355 * abs(sx_f)
        fwhm_y = 2.355 * abs(sy_f)
        return fwhm_x, fwhm_y, float(theta_f)
    except Exception:
        # Fall back to moment-based estimate
        return 2.355 * sx, 2.355 * sy, 0.0


# ─────────────────────────────────────────────────────────────────────────────
# ePSF → native pixel-scale kernel
# ─────────────────────────────────────────────────────────────────────────────

def _epsf_to_kernel(epsf_model, oversampling, output_size=None):
    """
    Extract the ePSF array from a photutils EPSFModel, downsample to native
    pixel scale, ensure odd size, normalise to unit sum.

    Parameters
    ----------
    epsf_model   : photutils EPSFModel
    oversampling : int used when building the ePSF
    output_size  : int or None — desired output kernel side length (odd).
                   If None, uses the full downsampled array.

    Returns
    -------
    kernel : 2-D ndarray, shape (output_size, output_size), sum=1
    """
    data = epsf_model.data.copy()
    data = np.clip(data, 0, None)

    # Downsample from oversampled to native pixel scale
    if oversampling > 1:
        data = zoom(data, 1.0 / oversampling, order=3, prefilter=True)
        data = np.clip(data, 0, None)

    # Trim or pad to requested output_size
    if output_size is not None:
        ny, nx = data.shape
        if output_size % 2 == 0:
            output_size += 1
        # Centre-crop
        cy, cx = ny // 2, nx // 2
        h = output_size // 2
        y1, y2 = max(0, cy - h), min(ny, cy + h + 1)
        x1, x2 = max(0, cx - h), min(nx, cx + h + 1)
        data = data[y1:y2, x1:x2]
        # Pad if crop was too small
        if data.shape[0] < output_size or data.shape[1] < output_size:
            pad_y = max(0, output_size - data.shape[0])
            pad_x = max(0, output_size - data.shape[1])
            data = np.pad(data, ((pad_y // 2, pad_y - pad_y // 2),
                                  (pad_x // 2, pad_x - pad_x // 2)))

    # Normalise
    s = data.sum()
    if s > 0:
        data /= s

    return data


# ─────────────────────────────────────────────────────────────────────────────
# PSF quality metric: FWHM scatter + peak-normalised residual RMS
# ─────────────────────────────────────────────────────────────────────────────

def _compute_psf_quality(native_kernel, star_cutouts,
                         saturation=55000, gain=1.77, rdnoise=4.0):
    """
    Compute two PSF quality metrics from the fitted star cutouts.

    Returns a dict with:
        fwhm_scatter : float
            Standard deviation of per-star FWHM values (px).  Low scatter
            means the PSF is consistent across the field (good model).
        residual_rms_pct : float
            Mean peak-normalised RMS of (data − model) in the PSF core,
            expressed as a percentage.  Evaluated only for stars in a
            moderate brightness range to avoid saturation and readout noise
            domination.  Independent of sub-pixel centering because model
            is aligned by shifting to match each star's measured centroid.

    Parameters
    ----------
    native_kernel : 2-D ndarray (sum≈1, native pixel scale)
    star_cutouts  : EPSFStars / list of EPSFStar returned by EPSFBuilder
    saturation    : float — ADU ceiling; stars above this are excluded
    gain, rdnoise : detector parameters (e-/ADU, e-)
    """
    if native_kernel is None or star_cutouts is None:
        return dict(fwhm_scatter=np.nan, residual_rms_pct=np.nan)

    kny, knx = native_kernel.shape
    kernel_norm = native_kernel / native_kernel.max()

    per_star_fwhm = []
    per_star_rms  = []

    for star in star_cutouts:
        data = np.array(star.data, dtype=float)
        ny, nx = data.shape

        # Background: median of corner pixels
        margin = max(2, ny // 6)
        corners = np.concatenate([
            data[:margin, :margin].ravel(), data[:margin, -margin:].ravel(),
            data[-margin:, :margin].ravel(), data[-margin:, -margin:].ravel(),
        ])
        bkg = float(np.median(corners))
        data = data - bkg
        peak = float(data.max())

        if peak < 200:            # too faint — Gaussian fit unreliable
            continue

        # ── per-star Gaussian FWHM (unsaturated stars only) ─────────────────
        if peak <= saturation:
            try:
                fwhm_xi, fwhm_yi, _ = _fit_psf_fwhm(data)
                fwhm_i = float(np.sqrt(fwhm_xi * fwhm_yi))
                if 0.5 < fwhm_i < 20:  # sanity range
                    per_star_fwhm.append(fwhm_i)
            except Exception:
                pass

        # ── peak-normalised residual (unsaturated, moderate-brightness only) ─
        if peak > saturation:
            continue

        data_norm = data / peak

        # Shift kernel to match star's measured centroid rather than assuming
        # the star is perfectly centred in the cutout.
        peak_iy, peak_ix = np.unravel_index(np.argmax(data), data.shape)
        shift_y = ny // 2 - peak_iy
        shift_x = nx // 2 - peak_ix
        if (ny, nx) == (kny, knx):
            model_norm = kernel_norm.copy()
        else:
            model_norm = zoom(kernel_norm, (ny / kny, nx / knx),
                              order=3, prefilter=True)
            model_norm = np.clip(model_norm, 0, None)
            if model_norm.max() > 0:
                model_norm /= model_norm.max()

        # Apply integer-pixel shift
        if shift_y != 0 or shift_x != 0:
            model_norm = np.roll(model_norm, (shift_y, shift_x), axis=(0, 1))

        # Central region mask (≤ 1.5 FWHM radius)
        cy, cx = ny / 2.0, nx / 2.0
        yy, xx = np.mgrid[:ny, :nx]
        r = np.hypot(xx - cx, yy - cy)
        mask = r <= max(kny // 3, 4)
        if mask.sum() < 9:
            continue

        rms = float(np.std(data_norm[mask] - model_norm[mask]))
        per_star_rms.append(rms * 100.0)   # convert to %

    # Sigma-clip FWHM list (1.5×IQR) to reject catastrophic outliers before scatter
    if len(per_star_fwhm) >= 4:
        fwhm_arr = np.array(per_star_fwhm)
        q1, q3 = np.percentile(fwhm_arr, [25, 75])
        iqr = q3 - q1
        lo, hi = q1 - 2.5 * iqr, q3 + 2.5 * iqr
        fwhm_arr = fwhm_arr[(fwhm_arr >= lo) & (fwhm_arr <= hi)]
        fwhm_scatter = float(np.std(fwhm_arr)) if len(fwhm_arr) >= 3 else np.nan
    elif len(per_star_fwhm) >= 3:
        fwhm_scatter = float(np.std(per_star_fwhm))
    else:
        fwhm_scatter = np.nan

    residual_rms_pct = float(np.mean(per_star_rms)) if per_star_rms else np.nan

    return dict(fwhm_scatter=fwhm_scatter, residual_rms_pct=residual_rms_pct)


# ─────────────────────────────────────────────────────────────────────────────
# Main API
# ─────────────────────────────────────────────────────────────────────────────

def build_psf(
    image_data,
    fwhm_guess=5.0,
    threshold_sigma=5.0,
    saturation=None,
    min_snr=10.0,
    oversampling=4,
    max_stars=40,
    min_stars=5,
    isolation_factor=2.5,
    kernel_size=None,
    bkg_box=64,
    epsf_iters=3,
    epsf_smoothing_kernel='quartic',
):
    """
    Build an empirical PSF from an astronomical image using photutils EPSFBuilder.

    Parameters
    ----------
    image_data : 2-D ndarray
    fwhm_guess : float
        Expected FWHM in pixels.
    threshold_sigma : float
        Detection threshold in units of background RMS.
    saturation : float or None
        ADU value above which pixels are saturated.
    min_snr : float
        Minimum peak SNR a star must have to be included in ePSF construction.
        Stars below this are noise spikes that corrupt the stacked PSF; they
        are silently dropped.  The SNR is estimated as peak / local_rms.
        Set to 0 to disable (not recommended).  Default 10.
    oversampling : int
        EPSFBuilder oversampling factor.  Higher = more accurate sub-pixel
        centering but slower. Typically 4 for SEDM pixel scale.
    max_stars : int
        Maximum number of (brightest isolated) stars to use.
    min_stars : int
        Minimum stars required; raises ValueError if not met.
    isolation_factor : float
        Stars with neighbours within isolation_factor * fwhm_guess are excluded.
    kernel_size : int or None
        Side length (odd, pixels) of the returned native-scale kernel.
        Defaults to round(10 * fwhm_guess) rounded to next odd integer
        (minimum 51). Using 10×FWHM captures the full stellar wing profile
        and matches PSFEx's default 59×59 kernel size at typical SEDM pixel scale.
    bkg_box : int
        Box size for 2-D background estimation.
    epsf_iters : int
        Number of EPSFBuilder fitting iterations.
    epsf_smoothing_kernel : str or 2-D array
        Smoothing kernel applied during ePSF construction. 'quartic' (default)
        is the photutils recommendation.

    Returns
    -------
    dict — see module docstring for full key list.
    """
    data = np.asarray(image_data, dtype=float)
    ny, nx = data.shape

    # ── default kernel output size ────────────────────────────────────────────
    # Use 10×FWHM to capture stellar wings beyond the PSF core.  6×FWHM
    # (old default) only extended to r≈13 px for typical FWHM≈4.5 px, causing
    # the radial profile to turn up relative to PSFEx past ~6 px because the
    # truncated wings were renormalised.  10×FWHM matches PSFEx's 59×59 kernel
    # for the typical SEDM pixel scale.
    if kernel_size is None:
        kernel_size = int(np.ceil(fwhm_guess * 10))
        if kernel_size % 2 == 0:
            kernel_size += 1
        kernel_size = max(kernel_size, 51)

    # ── 1. Background subtraction ─────────────────────────────────────────────
    bkg_2d, rms_2d = _estimate_background(data, box_size=bkg_box)
    data_sub = data - bkg_2d

    # ── 2. Detect sources ─────────────────────────────────────────────────────
    sources = _detect_stars(
        data_sub, rms_2d,
        fwhm_guess=fwhm_guess,
        threshold_sigma=threshold_sigma,
        saturation=saturation,
        edge_margin=max(20, kernel_size // 2 + 5),
    )

    if len(sources) == 0:
        raise ValueError(
            f"No sources detected above {threshold_sigma}σ. "
            "Try lowering threshold_sigma or adjusting fwhm_guess."
        )

    # ── 3. Isolation filter ───────────────────────────────────────────────────
    sources_iso = _filter_isolated(sources, fwhm=fwhm_guess,
                                   isolation_factor=isolation_factor)

    if len(sources_iso) < min_stars:
        # Relax isolation if we don't have enough
        sources_iso = _filter_isolated(sources, fwhm=fwhm_guess,
                                       isolation_factor=1.8)
        if len(sources_iso) < min_stars:
            raise ValueError(
                f"Only {len(sources_iso)} isolated stars found (need {min_stars}). "
                "Try lowering min_stars or threshold_sigma."
            )

    # ── 3b. SNR floor: drop near-noise detections ─────────────────────────────
    # Stars with peak < min_snr × local_rms have no meaningful PSF signal;
    # including them corrupts the stacked ePSF with noise-dominated cutouts.
    if min_snr > 0:
        local_rms = float(np.median(rms_2d))
        snr_mask = sources_iso['peak'] >= min_snr * local_rms
        sources_snr = sources_iso[snr_mask]
        if len(sources_snr) >= min_stars:
            sources_iso = sources_snr
        # else: keep all stars (better noisy PSF than no PSF)

    # Sort by flux, take the brightest max_stars
    sources_iso.sort('flux', reverse=True)
    sources_iso = sources_iso[:max_stars]

    # ── 4. Extract star cutouts ───────────────────────────────────────────────
    # EPSFBuilder wants cutout size in native pixels
    cutout_size = kernel_size  # odd integer

    stars_tbl = Table()
    stars_tbl['x'] = sources_iso['xcentroid']
    stars_tbl['y'] = sources_iso['ycentroid']

    nddata = NDData(data=data_sub)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        stars = extract_stars(nddata, stars_tbl, size=cutout_size)

    if len(stars) < min_stars:
        raise ValueError(
            f"extract_stars returned only {len(stars)} cutouts (need {min_stars})."
        )

    # ── 5. Build ePSF ─────────────────────────────────────────────────────────
    # photutils 1.x uses EPSFBuilder()(stars), 2.x uses EPSFBuilder().build_epsf(stars)
    builder_kwargs = dict(
        oversampling=oversampling,
        maxiters=epsf_iters,
        progress_bar=False,
        smoothing_kernel=epsf_smoothing_kernel,
    )
    # recentering_maxiters exists in 1.x but not 2.x
    try:
        import inspect as _inspect
        if 'recentering_maxiters' in _inspect.signature(EPSFBuilder.__init__).parameters:
            builder_kwargs['recentering_maxiters'] = 10
    except Exception:
        pass

    builder = EPSFBuilder(**builder_kwargs)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        # Try photutils 2.x API first, fall back to 1.x
        if hasattr(builder, 'build_epsf'):
            epsf, fitted_stars = builder.build_epsf(stars)
        else:
            epsf, fitted_stars = builder(stars)

    # ── 6. Convert to native pixel-scale kernel ───────────────────────────────
    kernel = _epsf_to_kernel(epsf, oversampling=oversampling,
                             output_size=kernel_size)

    # ── 7. FWHM from Gaussian fit to kernel ───────────────────────────────────
    fwhm_x, fwhm_y, theta = _fit_psf_fwhm(kernel)
    fwhm = float(np.sqrt(fwhm_x * fwhm_y))
    elongation = (max(fwhm_x, fwhm_y) /
                  max(min(fwhm_x, fwhm_y), 1e-6))

    # ── 8. PSF quality metrics ────────────────────────────────────────────────
    quality = _compute_psf_quality(kernel, fitted_stars)

    return dict(
        kernel           = kernel,
        epsf             = epsf,
        fwhm             = fwhm,
        fwhm_x           = fwhm_x,
        fwhm_y           = fwhm_y,
        theta            = theta,
        elongation       = elongation,
        # Legacy key kept for back-compat; now holds FWHM scatter (px)
        chi2             = quality['fwhm_scatter'],
        fwhm_scatter     = quality['fwhm_scatter'],
        residual_rms_pct = quality['residual_rms_pct'],
        n_stars          = len(fitted_stars),
        stars_x          = np.array(sources_iso['xcentroid'][:len(fitted_stars)]),
        stars_y          = np.array(sources_iso['ycentroid'][:len(fitted_stars)]),
        oversampling     = oversampling,
    )


def build_psf_from_fits(filepath, ext=0, **kwargs):
    """
    Load a FITS file and build its PSF.

    Parameters
    ----------
    filepath : str or Path
    ext      : int or str — FITS extension to read
    **kwargs : passed to build_psf()

    Returns
    -------
    Same dict as build_psf().
    """
    with fits.open(filepath) as hdul:
        data = hdul[ext].data
        header = hdul[ext].header

    if data is None:
        raise ValueError(f"No image data in extension {ext} of {filepath}")

    # Auto-read saturation from header if not provided
    if 'saturation' not in kwargs:
        sat = header.get('SATURATE', header.get('SATURLEV', None))
        if sat is not None:
            kwargs['saturation'] = float(sat) * 0.95   # 5% margin

    return build_psf(data.astype(float), **kwargs)
