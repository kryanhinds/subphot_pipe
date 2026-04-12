"""
psf_measure.py — Stable empirical PSF measurement for astronomical images.

Requires: numpy, scipy, astropy

Usage
-----
from psf_measure import measure_psf, measure_psf_from_fits

# From a numpy array (raw ADU image)
result = measure_psf(data, fwhm_guess=5.0, saturation=60000)

# From a FITS file
result = measure_psf_from_fits('image.fits', fwhm_guess=5.0)

print(f"FWHM = {result['fwhm']:.2f} px, using {result['n_stars']} stars")

Returns dict with keys:
  'psf'        2D ndarray — median-stacked, flux-normalised PSF stamp
  'fwhm'       float      — circularised FWHM = sqrt(fwhm_x * fwhm_y)  [pixels]
  'fwhm_x'     float      — FWHM along x-axis  [pixels]
  'fwhm_y'     float      — FWHM along y-axis  [pixels]
  'theta'      float      — PSF position angle  [radians]
  'elongation' float      — major/minor FWHM ratio
  'n_stars'    int        — number of stars used
  'stars_x'    ndarray    — x centroids of PSF stars
  'stars_y'    ndarray    — y centroids of PSF stars
"""

import warnings
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from scipy.ndimage import gaussian_filter, maximum_filter
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import curve_fit


# ─────────────────────────────────────────────────────────────────────────────
# 2-D profile model
# ─────────────────────────────────────────────────────────────────────────────

def _gaussian_2d(xy, amplitude, x0, y0, sigma_x, sigma_y, theta, offset):
    """
    Tilted 2-D Gaussian, returns flattened array.

    exponent = -0.5 * (a*dx^2 + b*dx*dy + c*dy^2)
    where the rotation is encoded in a, b, c.
    """
    x, y = xy
    dx, dy = x - x0, y - y0
    ct, st = np.cos(theta), np.sin(theta)
    a = (ct / sigma_x) ** 2 + (st / sigma_y) ** 2
    b = np.sin(2 * theta) * (1 / sigma_x ** 2 - 1 / sigma_y ** 2)
    c = (st / sigma_x) ** 2 + (ct / sigma_y) ** 2
    z = offset + amplitude * np.exp(-0.5 * (a * dx**2 + b * dx * dy + c * dy**2))
    return z.ravel()


# ─────────────────────────────────────────────────────────────────────────────
# Background estimation
# ─────────────────────────────────────────────────────────────────────────────

def _sigma_clip_median(arr, sigma=3.0, maxiters=5):
    """Iterative sigma-clip; return median and NMAD-based RMS."""
    arr = arr[np.isfinite(arr)].ravel()
    for _ in range(maxiters):
        med = np.median(arr)
        mad = np.median(np.abs(arr - med))
        rms = 1.4826 * mad          # NMAD → equivalent Gaussian sigma
        if rms == 0:
            break
        arr = arr[np.abs(arr - med) <= sigma * rms]
    med = float(np.median(arr)) if len(arr) else 0.0
    rms = 1.4826 * float(np.median(np.abs(arr - med))) if len(arr) else 1.0
    return med, rms


def _estimate_background(image, box_size=64, clip_sigma=3.0):
    """
    Estimate a smooth 2-D background by sigma-clipped median on a grid,
    then bicubic-spline interpolation.

    Returns
    -------
    bkg : 2-D array, same shape as image
    rms : float, robust background RMS
    """
    ny, nx = image.shape
    step = max(16, min(box_size, ny // 4, nx // 4))

    # Grid of box centres
    ys = np.arange(step // 2, ny, step)
    xs = np.arange(step // 2, nx, step)

    # Clamp so we always have at least one point per axis
    ys = np.unique(np.clip(ys, 0, ny - 1))
    xs = np.unique(np.clip(xs, 0, nx - 1))

    grid = np.zeros((len(ys), len(xs)))
    for iy, y0 in enumerate(ys):
        for ix, x0 in enumerate(xs):
            box = image[max(0, y0 - step // 2): y0 + step // 2 + 1,
                        max(0, x0 - step // 2): x0 + step // 2 + 1]
            grid[iy, ix], _ = _sigma_clip_median(box.ravel(), sigma=clip_sigma)

    # Interpolate to full resolution
    if len(ys) >= 4 and len(xs) >= 4:
        try:
            spline = RectBivariateSpline(ys, xs, grid, kx=min(3, len(ys)-1),
                                         ky=min(3, len(xs)-1))
            bkg = spline(np.arange(ny), np.arange(nx))
        except Exception:
            bkg = np.full((ny, nx), np.median(grid))
    else:
        bkg = np.full((ny, nx), np.median(grid))

    # Global RMS on residual
    _, rms = _sigma_clip_median((image - bkg).ravel(), sigma=clip_sigma)
    rms = max(rms, 1e-6)   # guard against zero-noise test images
    return bkg, rms


# ─────────────────────────────────────────────────────────────────────────────
# Source detection
# ─────────────────────────────────────────────────────────────────────────────

def _find_peaks(image_sub, threshold, min_sep):
    """
    Find local maxima above threshold with a minimum separation.

    Parameters
    ----------
    image_sub  : background-subtracted 2-D array
    threshold  : minimum peak value
    min_sep    : minimum distance between peaks [pixels]

    Returns
    -------
    List of (x, y, peak_value) sorted brightest-first.
    """
    neigh = max(3, int(2 * min_sep + 1))
    if neigh % 2 == 0:
        neigh += 1

    local_max = maximum_filter(image_sub, size=neigh)
    peak_mask = (image_sub == local_max) & (image_sub > threshold)

    ys, xs = np.where(peak_mask)
    if len(xs) == 0:
        return []

    vals = image_sub[ys, xs]
    order = np.argsort(vals)[::-1]
    xs, ys, vals = xs[order], ys[order], vals[order]

    # Greedy deblend: drop peaks within min_sep of a brighter accepted peak
    keep = np.ones(len(xs), dtype=bool)
    for i in range(len(xs)):
        if not keep[i]:
            continue
        for j in range(i + 1, len(xs)):
            if not keep[j]:
                continue
            if np.hypot(xs[i] - xs[j], ys[i] - ys[j]) < min_sep:
                keep[j] = False

    return list(zip(xs[keep].tolist(), ys[keep].tolist(), vals[keep].tolist()))


def _refine_centroid(image_sub, x0, y0, radius):
    """Intensity-weighted centroid within a circular aperture."""
    ny, nx = image_sub.shape
    r = int(np.ceil(radius))
    x1, x2 = max(0, int(x0) - r), min(nx, int(x0) + r + 1)
    y1, y2 = max(0, int(y0) - r), min(ny, int(y0) + r + 1)

    sub = image_sub[y1:y2, x1:x2].copy()
    sub = np.clip(sub - np.min(sub), 0, None)

    total = sub.sum()
    if total == 0:
        return x0, y0

    ys_g, xs_g = np.mgrid[y1:y2, x1:x2]
    return float((xs_g * sub).sum() / total), float((ys_g * sub).sum() / total)


# ─────────────────────────────────────────────────────────────────────────────
# Per-star Gaussian fitting
# ─────────────────────────────────────────────────────────────────────────────

def _fit_star(cutout):
    """
    Fit a tilted 2-D Gaussian to a background-subtracted star cutout.

    Returns (sigma_x, sigma_y, theta, success).
    Sizes are in pixels; sigma → FWHM via 2.355 * sigma.
    """
    ny, nx = cutout.shape
    y_g, x_g = np.mgrid[:ny, :nx]

    # Robust background: 10th-percentile of cutout
    bkg = float(np.percentile(cutout, 10))
    data = np.clip(cutout - bkg, 0, None)

    peak = float(data.max())
    if peak <= 0:
        return None, None, None, False

    total = data.sum()
    cx = float((x_g * data).sum() / total) if total > 0 else nx / 2.0
    cy = float((y_g * data).sum() / total) if total > 0 else ny / 2.0
    dx = x_g - cx
    dy = y_g - cy
    sx = float(np.sqrt(max(0.25, (dx**2 * data).sum() / total)))
    sy = float(np.sqrt(max(0.25, (dy**2 * data).sum() / total)))

    p0 = [peak, cx, cy, sx, sy, 0.0, bkg]
    lim = min(nx, ny)
    bounds = (
        [0,   -1,  -1,  0.3, 0.3, -np.pi / 2, -np.inf],
        [np.inf, nx+1, ny+1, lim, lim,  np.pi / 2,  np.inf],
    )

    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            popt, _ = curve_fit(
                _gaussian_2d, (x_g, y_g), cutout.ravel(),
                p0=p0, bounds=bounds, maxfev=3000,
                ftol=1e-6, xtol=1e-6,
            )
        amp, x0, y0, sigma_x, sigma_y, theta, _ = popt

        # Sanity checks
        if amp <= 0:
            return None, None, None, False
        if min(sigma_x, sigma_y) < 0.3 or max(sigma_x, sigma_y) > lim * 0.9:
            return None, None, None, False
        if abs(x0 - nx / 2) > nx / 2 or abs(y0 - ny / 2) > ny / 2:
            return None, None, None, False

        return abs(sigma_x), abs(sigma_y), float(theta), True

    except Exception:
        return None, None, None, False


# ─────────────────────────────────────────────────────────────────────────────
# Main API
# ─────────────────────────────────────────────────────────────────────────────

def measure_psf(
    image_data,
    fwhm_guess=5.0,
    threshold_sigma=5.0,
    saturation=None,
    edge_margin=None,
    box_size=None,
    min_stars=3,
    max_stars=50,
    bkg_box=64,
):
    """
    Measure the PSF of an astronomical image.

    The algorithm:
    1. Estimate a smooth 2-D background via sigma-clipped grid medians.
    2. Detect point sources as local maxima above threshold_sigma × rms.
    3. Filter candidates: edge proximity, saturation, compactness (star/galaxy).
    4. Fit a tilted 2-D Gaussian to each star cutout.
    5. Sigma-clip FWHM measurements; report robust median.
    6. Build an empirical PSF by median-stacking normalised star stamps.

    Parameters
    ----------
    image_data : 2-D ndarray
        Raw image pixel values (any units).
    fwhm_guess : float
        Expected stellar FWHM in pixels.  Used to set the detection smoothing
        kernel and cutout size.  A ±2× error is tolerated.
    threshold_sigma : float
        Detection threshold in multiples of background RMS.
    saturation : float or None
        ADU value above which a pixel is saturated.  Stars whose cutout peak
        exceeds this are excluded.
    edge_margin : int or None
        Exclusion zone at image edges [pixels].  Default = box_size // 2 + 3.
    box_size : int or None
        Side length of the PSF stamp [pixels, odd].  Default = 8 × fwhm_guess,
        rounded to the next odd integer, minimum 21.
    min_stars : int
        Minimum usable stars; raises ValueError if not met.
    max_stars : int
        Cap on the number of (brightest) stars used.
    bkg_box : int
        Grid-box size for background estimation.

    Returns
    -------
    dict with keys:
        'psf'        : 2-D ndarray — flux-normalised empirical PSF stamp
        'fwhm'       : float — circularised FWHM = sqrt(fwhm_x × fwhm_y)  [px]
        'fwhm_x'     : float — FWHM along x  [px]
        'fwhm_y'     : float — FWHM along y  [px]
        'theta'      : float — position angle of PSF major axis  [radians]
        'elongation' : float — major / minor FWHM ratio  (1 = round)
        'n_stars'    : int   — number of stars used
        'stars_x'    : ndarray — x centroids of used stars
        'stars_y'    : ndarray — y centroids of used stars

    Raises
    ------
    ValueError if fewer than min_stars usable stars are found.
    """
    image_data = np.asarray(image_data, dtype=float)
    ny, nx = image_data.shape

    # ── Derived defaults ──────────────────────────────────────────────────────
    if box_size is None:
        box_size = int(np.ceil(fwhm_guess * 8))
        if box_size % 2 == 0:
            box_size += 1
        box_size = max(box_size, 21)
    elif box_size % 2 == 0:
        box_size += 1

    half = box_size // 2

    if edge_margin is None:
        edge_margin = half + 3

    sigma_detect = fwhm_guess / 2.355

    # ── 1. Background ─────────────────────────────────────────────────────────
    bkg, bkg_rms = _estimate_background(image_data, box_size=bkg_box)
    image_sub = image_data - bkg

    # ── 2. Detect peaks on a Gaussian-smoothed image ──────────────────────────
    smoothed = gaussian_filter(image_sub, sigma=sigma_detect)
    threshold = threshold_sigma * bkg_rms
    peaks = _find_peaks(smoothed, threshold, min_sep=fwhm_guess)

    if not peaks:
        raise ValueError(
            f"No sources detected above {threshold_sigma}σ "
            f"(threshold = {threshold:.1f} ADU). "
            "Try lowering threshold_sigma or providing a better fwhm_guess."
        )

    # ── 3. Filter candidates ──────────────────────────────────────────────────
    candidates = []
    inner_r = fwhm_guess / 2.0       # ~half-light radius for a star
    outer_r = fwhm_guess * 2.0       # encircles >99 % of star flux

    for (px, py, peak_smooth) in peaks:
        # Coarse edge check before centroid refinement
        if not (edge_margin <= px < nx - edge_margin and
                edge_margin <= py < ny - edge_margin):
            continue

        # Refine centroid
        cx, cy = _refine_centroid(image_sub, px, py, radius=fwhm_guess)

        # Strict edge check with box half-size
        if not (half < cx < nx - half and half < cy < ny - half):
            continue

        # Integer pixel for array slicing
        xi, yi = int(round(cx)), int(round(cy))
        x1, x2 = xi - half, xi + half + 1
        y1, y2 = yi - half, yi + half + 1

        if x1 < 0 or y1 < 0 or x2 > nx or y2 > ny:
            continue

        cutout_raw = image_data[y1:y2, x1:x2]
        cutout_sub = image_sub[y1:y2, x1:x2]

        # Saturation
        if saturation is not None and np.max(cutout_raw) > saturation:
            continue

        # SNR check on un-smoothed residual
        peak_val = float(np.max(cutout_sub))
        if peak_val < threshold:
            continue

        # Compactness (star / galaxy discriminator).
        # For a Gaussian PSF, ~50 % of flux is within r < fwhm/2.
        # We require the inner-to-outer flux ratio > 0.25 (quite loose).
        cy_c = box_size / 2.0
        cx_c = box_size / 2.0
        yy, xx = np.mgrid[:box_size, :box_size]
        rr = np.hypot(xx - cx_c, yy - cy_c)
        s_inner = float(cutout_sub[rr <= inner_r].sum())
        s_outer = float(cutout_sub[rr <= outer_r].sum())
        if s_outer <= 0 or (s_inner / s_outer) < 0.20:
            continue   # too extended → galaxy or artefact

        candidates.append(dict(
            x=cx, y=cy,
            peak=peak_val,
            cutout_sub=cutout_sub.copy(),
        ))

    if not candidates:
        raise ValueError(
            "No candidates survived quality filters (edge, saturation, compactness). "
            "Try lowering threshold_sigma or saturation."
        )

    # Sort brightest first, cap at max_stars
    candidates.sort(key=lambda s: s['peak'], reverse=True)
    candidates = candidates[:max_stars]

    # ── 4. Fit 2-D Gaussian to each candidate ─────────────────────────────────
    fitted = []
    psf_stamps = []

    fwhm_lo = fwhm_guess * 0.25   # floor: reject cosmic rays / hot pixels
    fwhm_hi = fwhm_guess * 4.0    # ceiling: reject very extended sources

    for cand in candidates:
        sx, sy, theta, ok = _fit_star(cand['cutout_sub'])
        if not ok:
            continue

        fwhm_x_i = 2.355 * sx
        fwhm_y_i = 2.355 * sy
        fwhm_i   = np.sqrt(fwhm_x_i * fwhm_y_i)

        if not (fwhm_lo <= fwhm_i <= fwhm_hi):
            continue

        # Normalised stamp for PSF building
        stamp = cand['cutout_sub'].copy()
        stamp_peak = float(stamp.max())
        if stamp_peak > 0:
            psf_stamps.append(stamp / stamp_peak)

        fitted.append(dict(
            x=cand['x'], y=cand['y'],
            fwhm_x=fwhm_x_i, fwhm_y=fwhm_y_i, theta=theta,
        ))

    if len(fitted) < min_stars:
        raise ValueError(
            f"Only {len(fitted)} stars passed the Gaussian fit quality checks "
            f"(need {min_stars}). "
            "Try lowering min_stars, adjusting fwhm_guess, or lowering threshold_sigma."
        )

    # ── 5. Sigma-clip FWHM across stars ───────────────────────────────────────
    def _sc(arr, sigma=2.5):
        med = np.median(arr)
        std = np.std(arr)
        return arr[np.abs(arr - med) <= sigma * std] if std > 0 else arr

    fwhm_x_arr = np.array([f['fwhm_x'] for f in fitted])
    fwhm_y_arr = np.array([f['fwhm_y'] for f in fitted])
    thetas     = np.array([f['theta']  for f in fitted])

    fwhm_x_final = float(np.median(_sc(fwhm_x_arr)))
    fwhm_y_final = float(np.median(_sc(fwhm_y_arr)))
    theta_final  = float(np.median(thetas))
    fwhm_final   = np.sqrt(fwhm_x_final * fwhm_y_final)
    elongation   = (max(fwhm_x_final, fwhm_y_final) /
                    max(min(fwhm_x_final, fwhm_y_final), 1e-6))

    # ── 6. Median-stack empirical PSF ──────────────────────────────────────────
    if psf_stamps:
        psf_array = np.median(psf_stamps, axis=0)
        psf_array = np.clip(psf_array, 0, None)
        s = psf_array.sum()
        if s > 0:
            psf_array /= s
    else:
        psf_array = None

    return dict(
        psf        = psf_array,
        fwhm       = fwhm_final,
        fwhm_x     = fwhm_x_final,
        fwhm_y     = fwhm_y_final,
        theta      = theta_final,
        elongation = elongation,
        n_stars    = len(fitted),
        stars_x    = np.array([f['x'] for f in fitted]),
        stars_y    = np.array([f['y'] for f in fitted]),
    )


def measure_psf_from_fits(filepath, ext=0, **kwargs):
    """
    Load a FITS file and measure the PSF.

    Parameters
    ----------
    filepath : str or Path
    ext      : int or str — FITS extension to read
    **kwargs : passed to measure_psf()

    Returns
    -------
    Same dict as measure_psf().
    """
    with fits.open(filepath) as hdul:
        data = hdul[ext].data
    if data is None:
        raise ValueError(f"No image data in extension {ext} of {filepath}")
    return measure_psf(data.astype(float), **kwargs)


# ─────────────────────────────────────────────────────────────────────────────
# Optional visualisation (requires matplotlib)
# ─────────────────────────────────────────────────────────────────────────────

def plot_psf_result(image_data, result, vmin_pct=1, vmax_pct=99):
    """
    Quick diagnostic plot: image with star positions + PSF stamp + radial profile.

    Parameters
    ----------
    image_data : 2-D ndarray
    result     : dict returned by measure_psf()
    """
    try:
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
    except ImportError:
        print("matplotlib not available; skipping plot.")
        return

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    # ── Panel 1: image with detected star positions ───────────────────────────
    ax = axes[0]
    vmin = np.percentile(image_data, vmin_pct)
    vmax = np.percentile(image_data, vmax_pct)
    ax.imshow(image_data, origin='lower', cmap='gray',
              vmin=vmin, vmax=vmax, aspect='equal', interpolation='nearest')
    ax.scatter(result['stars_x'], result['stars_y'],
               s=60, facecolors='none', edgecolors='lime', linewidths=1.2)
    ax.set_title(f"Detected stars (n={result['n_stars']})")
    ax.set_xlabel('x  [px]')
    ax.set_ylabel('y  [px]')

    # ── Panel 2: empirical PSF stamp ─────────────────────────────────────────
    ax = axes[1]
    psf = result['psf']
    if psf is not None:
        p = np.clip(psf, psf[psf > 0].min() if (psf > 0).any() else 1e-9, None)
        ax.imshow(p, origin='lower', cmap='hot',
                  norm=LogNorm(vmin=p.min(), vmax=p.max()),
                  interpolation='nearest')
        ax.set_title('Empirical PSF (log scale)')
        ax.set_xlabel('x  [px]')
        ax.set_ylabel('y  [px]')
    else:
        ax.text(0.5, 0.5, 'PSF unavailable', ha='center', va='center',
                transform=ax.transAxes)

    # ── Panel 3: azimuthally averaged radial profile ─────────────────────────
    ax = axes[2]
    if psf is not None:
        ny, nx = psf.shape
        yc, xc = ny / 2.0, nx / 2.0
        yy, xx = np.mgrid[:ny, :nx]
        rr = np.hypot(xx - xc, yy - yc).ravel()
        zz = psf.ravel()
        order = np.argsort(rr)
        rr, zz = rr[order], zz[order]
        ax.plot(rr, zz, '.', ms=2, alpha=0.4, color='steelblue')

        # Binned median profile
        bins = np.arange(0, rr.max(), 0.5)
        bin_med = [np.median(zz[(rr >= bins[i]) & (rr < bins[i+1])])
                   for i in range(len(bins)-1)]
        ax.plot(bins[:-1] + 0.25, bin_med, 'r-', lw=1.5, label='median')

        half_peak = 0.5 * max(bin_med) if bin_med else 0
        ax.axhline(half_peak, ls='--', color='gray', lw=0.8, label='half-max')
        ax.axvline(result['fwhm'] / 2, ls=':', color='orange', lw=1.2,
                   label=f"FWHM/2 = {result['fwhm']/2:.2f} px")
        ax.set_yscale('log')
        ax.set_xlabel('Radius  [px]')
        ax.set_ylabel('Normalised flux')
        ax.set_title(f"Radial profile  (FWHM = {result['fwhm']:.2f} px)")
        ax.legend(fontsize=8)

    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# Quick smoke-test
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    import sys

    if len(sys.argv) > 1:
        path = sys.argv[1]
        fwhm = float(sys.argv[2]) if len(sys.argv) > 2 else 5.0
        print(f"Measuring PSF of {path}  (fwhm_guess={fwhm} px)...")
        res = measure_psf_from_fits(path, fwhm_guess=fwhm)
    else:
        # Synthetic test: 30 Gaussians on Poisson noise background
        rng = np.random.default_rng(42)
        ny, nx = 512, 512
        true_fwhm = 4.3
        true_sigma = true_fwhm / 2.355
        img = rng.poisson(200.0, (ny, nx)).astype(float)

        yy, xx = np.mgrid[:ny, :nx]
        for _ in range(30):
            x0 = rng.integers(60, nx - 60)
            y0 = rng.integers(60, ny - 60)
            amp = rng.uniform(3000, 15000)
            img += amp * np.exp(-0.5 * ((xx - x0)**2 + (yy - y0)**2) / true_sigma**2)

        print(f"Synthetic test: true FWHM = {true_fwhm:.2f} px")
        res = measure_psf(img, fwhm_guess=true_fwhm, threshold_sigma=5.0)

    print(f"  FWHM       = {res['fwhm']:.3f} px")
    print(f"  FWHM_x     = {res['fwhm_x']:.3f} px")
    print(f"  FWHM_y     = {res['fwhm_y']:.3f} px")
    print(f"  theta      = {np.degrees(res['theta']):.1f} deg")
    print(f"  elongation = {res['elongation']:.3f}")
    print(f"  n_stars    = {res['n_stars']}")
