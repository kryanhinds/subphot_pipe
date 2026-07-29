"""
subphot_psf.py — guaranteed PSF measurement for subphot_pipe.

One routine, `measure_psf`, that ALWAYS returns a usable,
normalised PSF kernel plus diagnostics — it never raises for data reasons
and never returns garbage silently.  Architecture (DAOPHOT / AutoPhOT /
LSST synthesis):

  C.0  Robust frame FWHM + shape bootstrap from second moments of 5-sigma
       detections (header seeing is only a cross-check, never trusted).
  C.1  Star selection: half-light-radius vs magnitude stellar locus
       (PSFEx/LSST style), SNR window, saturation margin, isolation,
       elongation cut RELATIVE to the frame median (legitimately
       elongated PSFs keep their stars), optional catalog cross-match.
  C.2  Stage 1 — sigma-clipped shifted star stack; per-star centroids from
       elliptical-Moffat fits (analytic recentring cannot run away the way
       EPSFBuilder's does), Piff-style worst-star rejection loop.
  C.3  Stage 2 — elliptical Moffat fit to the stack (converges with >=1 star).
  C.4  Stage 3 — hybrid = Moffat + capped, apodised residual table.
  C.5  Mode ladder:
         n>=5 good stars & clean residuals -> 'hybrid'
         2<=n<5                            -> 'moffat'
         n<2                               -> 'analytic'  (frame moments)
         no detections at all              -> 'fallback'  (seeing hint)

Every return carries a `diagnostics` dict (mode, n_stars, fwhm, residual
RMS, per-star scatter, warnings) so a low-quality PSF flags the
measurement downstream instead of killing it.

Usage
-----
from subphot_psf import measure_psf, measure_psf_from_fits

res = measure_psf(data, saturation=60000, seeing_hint_px=4.0)
kernel = res['kernel']          # 2-D ndarray, odd-sized, sum = 1
print(res['mode'], res['fwhm'], res['diagnostics']['warnings'])
"""

import warnings
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clip, sigma_clipped_stats, SigmaClip
from scipy import ndimage
from scipy.optimize import curve_fit

try:
    from photutils.background import Background2D, SExtractorBackground
    from photutils.detection import DAOStarFinder
    _HAS_PHOTUTILS = True
except ImportError:
    _HAS_PHOTUTILS = False

FWHM_MIN_PX, FWHM_MAX_PX = 1.5, 15.0      # hard clamps for fallback modes
_MOFFAT_BETA_DEFAULT = 2.5


# ─────────────────────────────────────────────────────────────────────────────
# Elliptical Moffat model
# ─────────────────────────────────────────────────────────────────────────────

def _moffat_2d(xy, amp, x0, y0, ax, ay, theta, beta, offset):
    """Elliptical Moffat; returns flattened array."""
    x, y = xy
    dx, dy = x - x0, y - y0
    ct, st = np.cos(theta), np.sin(theta)
    xr = dx * ct + dy * st
    yr = -dx * st + dy * ct
    u = (xr / ax) ** 2 + (yr / ay) ** 2
    return (offset + amp * (1.0 + u) ** (-beta)).ravel()


def _moffat_fwhm(alpha, beta):
    return 2.0 * alpha * np.sqrt(2.0 ** (1.0 / beta) - 1.0)


def _alpha_from_fwhm(fwhm, beta):
    return fwhm / (2.0 * np.sqrt(2.0 ** (1.0 / beta) - 1.0))


def _fit_moffat(stamp, fwhm_guess, fix_beta=None):
    """
    Fit an elliptical Moffat to a background-subtracted stamp.

    Returns dict(amp, x0, y0, ax, ay, theta, beta, offset,
                 fwhm_x, fwhm_y, fwhm, elongation, rms, ok)
    """
    ny, nx = stamp.shape
    y_g, x_g = np.mgrid[:ny, :nx]
    peak = float(np.nanmax(stamp))
    if not np.isfinite(peak) or peak <= 0:
        return dict(ok=False)

    data = np.nan_to_num(stamp, nan=0.0)
    beta0 = fix_beta if fix_beta is not None else _MOFFAT_BETA_DEFAULT
    a0 = max(0.7, _alpha_from_fwhm(fwhm_guess, beta0))
    iy, ix = np.unravel_index(np.argmax(data), data.shape)

    p0 = [peak, float(ix), float(iy), a0, a0, 0.0, beta0, 0.0]
    lim = max(nx, ny)
    lo = [0.0, 0.0, 0.0, 0.4, 0.4, -np.pi / 2, 1.2, -np.inf]
    hi = [np.inf, nx - 1, ny - 1, lim, lim, np.pi / 2, 10.0, np.inf]
    if fix_beta is not None:
        lo[6], hi[6] = beta0 - 1e-6, beta0 + 1e-6
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            popt, _ = curve_fit(_moffat_2d, (x_g, y_g), data.ravel(),
                                p0=p0, bounds=(lo, hi), maxfev=8000,
                                ftol=1e-7, xtol=1e-7)
    except Exception:
        return dict(ok=False)

    amp, x0, y0, ax, ay, theta, beta, offset = popt
    if amp <= 0 or min(ax, ay) < 0.3:
        return dict(ok=False)
    fwhm_x = _moffat_fwhm(abs(ax), beta)
    fwhm_y = _moffat_fwhm(abs(ay), beta)
    model = _moffat_2d((x_g, y_g), *popt).reshape(ny, nx)
    rms = float(np.sqrt(np.mean((data - model) ** 2))) / max(peak, 1e-9)
    return dict(amp=float(amp), x0=float(x0), y0=float(y0),
                ax=float(abs(ax)), ay=float(abs(ay)), theta=float(theta),
                beta=float(beta), offset=float(offset),
                fwhm_x=float(fwhm_x), fwhm_y=float(fwhm_y),
                fwhm=float(np.sqrt(fwhm_x * fwhm_y)),
                elongation=float(max(fwhm_x, fwhm_y) /
                                 max(min(fwhm_x, fwhm_y), 1e-6)),
                rms=rms, ok=True)


def _moffat_kernel(size, fwhm_x, fwhm_y, theta, beta):
    """Evaluate a unit-sum elliptical Moffat on an odd size x size grid."""
    if size % 2 == 0:
        size += 1
    c = size // 2
    y_g, x_g = np.mgrid[:size, :size]
    ax = _alpha_from_fwhm(max(fwhm_x, 0.8), beta)
    ay = _alpha_from_fwhm(max(fwhm_y, 0.8), beta)
    k = _moffat_2d((x_g, y_g), 1.0, c, c, ax, ay, theta, beta, 0.0)
    k = k.reshape(size, size)
    k = np.clip(k, 0, None)
    return k / k.sum()


# ─────────────────────────────────────────────────────────────────────────────
# Background + detection
# ─────────────────────────────────────────────────────────────────────────────

def _estimate_background(data, box_size=64):
    """2-D background + rms; photutils first, global sigma-clip fallback."""
    ny, nx = data.shape
    box = int(max(32, min(box_size, ny // 6 or 32, nx // 6 or 32)))
    if _HAS_PHOTUTILS:
        try:
            bkg = Background2D(data, box_size=(box, box), filter_size=(3, 3),
                               sigma_clip=SigmaClip(sigma=3.0, maxiters=5),
                               bkg_estimator=SExtractorBackground())
            return bkg.background, bkg.background_rms
        except Exception:
            pass
    _, med, rms = sigma_clipped_stats(data, sigma=3.0, maxiters=5)
    return (np.full_like(data, med, dtype=float),
            np.full_like(data, max(rms, 1e-6), dtype=float))


def _detect_sources(snr_img, threshold_sigma=5.0, fwhm_guesses=(3.0, 7.0),
                    max_sources=3000):
    """
    DAOStarFinder on a per-pixel SNR image (handles spatially varying
    noise and zero-padded borders from image alignment) at two FWHM
    guesses, merged.  Roundness/sharpness gates are opened wide — shape
    filtering happens later, relative to the frame's own PSF, so
    elongated-PSF frames are not starved of stars.
    Returns arrays (x, y, peak_snr), brightest-first, capped.
    """
    if not _HAS_PHOTUTILS:
        return (np.array([]),) * 3
    xs, ys, snrs = [], [], []
    for fw in fwhm_guesses:
        try:
            finder = DAOStarFinder(threshold=threshold_sigma,
                                   fwhm=fw, roundlo=-2.0, roundhi=2.0,
                                   sharplo=0.05, sharphi=2.0)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                tbl = finder(snr_img)
        except Exception:
            tbl = None
        if tbl is None or len(tbl) == 0:
            continue
        xs.extend(np.asarray(tbl['xcentroid'], float))
        ys.extend(np.asarray(tbl['ycentroid'], float))
        snrs.extend(np.asarray(tbl['peak'], float))
    if len(xs) == 0:
        return (np.array([]),) * 3
    x, y, sn = np.array(xs), np.array(ys), np.array(snrs)
    order = np.argsort(sn)[::-1][: 2 * max_sources]
    x, y, sn = x[order], y[order], sn[order]
    # dedupe detections from the two runs (within 2 px, keep brighter)
    keep = np.ones(len(x), bool)
    for i in range(len(x)):
        if not keep[i]:
            continue
        d = np.hypot(x - x[i], y - y[i])
        dup = (d < 2.0) & keep
        dup[i] = False
        keep[dup] = False
    x, y, sn = x[keep], y[keep], sn[keep]
    return x[:max_sources], y[:max_sources], sn[:max_sources]


def _neighbor_ratio(data_sub, x, y):
    """
    Second-brightest 4-neighbour / centre pixel.  A real PSF (FWHM >~ 1.8 px)
    puts >~0.35 of its peak in adjacent pixels; a cosmic-ray hit does not.
    Scale-free CR discriminator.
    """
    ny, nx = data_sub.shape
    xi, yi = int(round(x)), int(round(y))
    if not (1 <= xi < nx - 1 and 1 <= yi < ny - 1):
        return 0.0
    st = data_sub[yi - 1: yi + 2, xi - 1: xi + 2]
    iy, ix = np.unravel_index(np.argmax(st), st.shape)
    yi, xi = yi + iy - 1, xi + ix - 1
    if not (1 <= xi < nx - 1 and 1 <= yi < ny - 1):
        return 0.0
    c = float(data_sub[yi, xi])
    if c <= 0:
        return 0.0
    nb = sorted([data_sub[yi - 1, xi], data_sub[yi + 1, xi],
                 data_sub[yi, xi - 1], data_sub[yi, xi + 1]])
    return float(nb[2]) / c            # 2nd-brightest: robust to one hot column


def _mode_centers(vals, bin_w=0.4):
    """Ascending centres of the dense modes of a 1-D size distribution."""
    v = np.asarray(vals, float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return []
    if v.size < 8:
        return [float(np.median(v))]
    bins = np.arange(v.min(), v.max() + 2 * bin_w, bin_w)
    cnt, edges = np.histogram(v, bins=bins)
    cnt_s = np.convolve(cnt, [1, 2, 1], mode='same') if len(cnt) >= 3 else cnt
    centers = 0.5 * (edges[:-1] + edges[1:])
    thresh = max(4.0, 0.2 * cnt_s.max())
    out = []
    for i in range(len(cnt_s)):
        left = cnt_s[i - 1] if i > 0 else -1
        right = cnt_s[i + 1] if i < len(cnt_s) - 1 else -1
        if cnt_s[i] >= thresh and cnt_s[i] >= left and cnt_s[i] >= right:
            out.append(float(centers[i]))
    if not out:
        out = [float(centers[int(np.argmax(cnt_s))])]
    merged = []
    for c in out:
        if merged and c - merged[-1] < 2 * bin_w:
            continue
        merged.append(c)
    return merged


def _smallest_mode(vals, bin_w=0.4):
    """
    Centre + MAD of the LOWEST dense mode of a 1-D size distribution.
    Stars form the smallest-size population (nothing real is sharper than
    the PSF); in deep frames galaxies outnumber stars and drag the median
    up, so the seeing/locus estimate must be the low mode, not the median.
    """
    v = np.asarray(vals, float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return np.nan, np.nan
    med0 = float(np.median(v))
    mad0 = float(1.4826 * np.median(np.abs(v - med0)))
    if v.size < 8:
        return med0, mad0
    bins = np.arange(v.min(), v.max() + 2 * bin_w, bin_w)
    cnt, edges = np.histogram(v, bins=bins)
    cnt_s = np.convolve(cnt, [1, 2, 1], mode='same') if len(cnt) >= 3 else cnt
    centers = 0.5 * (edges[:-1] + edges[1:])
    thresh = max(4.0, 0.25 * cnt_s.max())
    mode_c = None
    for i in range(len(cnt_s)):
        left = cnt_s[i - 1] if i > 0 else -1
        right = cnt_s[i + 1] if i < len(cnt_s) - 1 else -1
        if cnt_s[i] >= thresh and cnt_s[i] >= left and cnt_s[i] >= right:
            mode_c = centers[i]                # first = smallest-value mode
            break
    if mode_c is None:
        mode_c = centers[int(np.argmax(cnt_s))]
    sel = v[np.abs(v - mode_c) <= max(2 * bin_w, 0.25 * abs(mode_c))]
    if sel.size == 0:
        return med0, mad0
    med = float(np.median(sel))
    return med, float(1.4826 * np.median(np.abs(sel - med)))


def _is_flattop(data, x, y, fwhm):
    """
    Unit-free saturation test: a saturated star shows a plateau of
    near-peak-equal pixels much larger than a real PSF core.  Works on
    RAW pixel values whatever the image units (raw ADU, gain-scaled,
    stacked floats), unlike any absolute ADU ceiling.
    """
    ny, nx = data.shape
    r = int(max(4, np.ceil(1.5 * fwhm)))
    xi, yi = int(round(x)), int(round(y))
    if not (r <= xi < nx - r and r <= yi < ny - r):
        return False
    st = data[yi - r: yi + r + 1, xi - r: xi + r + 1]
    peak = float(np.nanmax(st))
    if not np.isfinite(peak) or peak <= 0:
        return False
    lab, _ = ndimage.label(st >= 0.98 * peak)
    iy, ix = np.unravel_index(np.nanargmax(st), st.shape)
    n_flat = int(np.sum(lab == lab[iy, ix]))
    # a Moffat/Gaussian core has ~0.023*FWHM^2 px above 0.98*peak; allow 3x
    return n_flat > max(4.0, 0.07 * fwhm ** 2)


def _fwhm_halfmax(data_sub, x, y, rad=13):
    """
    FWHM from the area of the contiguous region above half the peak —
    immune to cosmic rays reading as bright 'stars' (their area is ~1 px)
    and to Moffat-wing inflation that biases moment-based FWHMs.
    """
    ny, nx = data_sub.shape
    xi, yi = int(round(x)), int(round(y))
    r = int(rad)
    if not (r <= xi < nx - r and r <= yi < ny - r):
        return np.nan
    st = data_sub[yi - r: yi + r + 1, xi - r: xi + r + 1]
    core = st[r - 1: r + 2, r - 1: r + 2]
    peak = float(np.nanmax(core))
    if not np.isfinite(peak) or peak <= 0:
        return np.nan
    mask = st >= 0.5 * peak
    lab, _ = ndimage.label(mask)
    lab_c = lab[r, r]
    if lab_c == 0:
        return np.nan
    area = float(np.sum(lab == lab_c))
    return 2.0 * np.sqrt(area / np.pi)


# ─────────────────────────────────────────────────────────────────────────────
# Per-source moments / half-light radius
# ─────────────────────────────────────────────────────────────────────────────

def _source_moments(data_sub, x, y, rad):
    """
    Annulus-background-subtracted second moments + half-flux radius in a
    circular window of radius `rad` around (x, y).
    Returns dict(fwhm, elongation, theta, rh, flux, snr_denom_area) or None.
    """
    ny, nx = data_sub.shape
    r = int(np.ceil(rad))
    xi, yi = int(round(x)), int(round(y))
    if not (r <= xi < nx - r and r <= yi < ny - r):
        return None
    sub = data_sub[yi - r: yi + r + 1, xi - r: xi + r + 1].astype(float)
    yy, xx = np.mgrid[:sub.shape[0], :sub.shape[1]]
    cy0, cx0 = r + (y - yi), r + (x - xi)
    rr = np.hypot(xx - cx0, yy - cy0)
    ann = sub[(rr > 0.75 * r) & (rr <= r)]
    if ann.size >= 8:
        local = float(np.median(ann))
        sub = sub - local
    win = rr <= 0.7 * r
    w = np.clip(sub, 0, None) * win
    tot = w.sum()
    if tot <= 0:
        return None
    cx = float((xx * w).sum() / tot)
    cy = float((yy * w).sum() / tot)
    dx, dy = xx - cx, yy - cy
    ixx = float((dx ** 2 * w).sum() / tot)
    iyy = float((dy ** 2 * w).sum() / tot)
    ixy = float((dx * dy * w).sum() / tot)
    t = ixx + iyy
    det = ixx * iyy - ixy ** 2
    if t <= 0 or det <= 0:
        return None
    lam1 = 0.5 * t + np.sqrt(max(0.25 * t ** 2 - det, 0))
    lam2 = 0.5 * t - np.sqrt(max(0.25 * t ** 2 - det, 0))
    if lam2 <= 0:
        return None
    fwhm = 2.355 * np.sqrt(0.5 * t)
    elong = np.sqrt(lam1 / lam2)
    theta = 0.5 * np.arctan2(2 * ixy, ixx - iyy)
    # half-flux radius from cumulative profile
    order = np.argsort(rr[win])
    csum = np.cumsum(np.clip(sub, 0, None)[win].ravel()[order])
    if csum[-1] <= 0:
        return None
    rh = float(np.sort(rr[win].ravel())[np.searchsorted(csum, 0.5 * csum[-1])])
    return dict(fwhm=float(fwhm), elongation=float(elong), theta=float(theta),
                rh=rh, flux=float(tot))


# ─────────────────────────────────────────────────────────────────────────────
# Kernel post-processing
# ─────────────────────────────────────────────────────────────────────────────

def _apodize(kernel, fwhm):
    """Cosine-taper the kernel to zero beyond ~3.5 x FWHM (ZTF-style
    wing regularisation) and clip tiny negatives."""
    ny, nx = kernel.shape
    c = (np.array(kernel.shape) - 1) / 2.0
    yy, xx = np.mgrid[:ny, :nx]
    rr = np.hypot(xx - c[1], yy - c[0])
    r0 = 3.0 * fwhm
    r1 = min(max(3.5 * fwhm, r0 + 2.0), 0.5 * min(nx, ny))
    taper = np.ones_like(kernel)
    m = rr > r0
    taper[m] = 0.5 * (1 + np.cos(np.pi * np.clip((rr[m] - r0) / max(r1 - r0, 1e-6), 0, 1)))
    out = kernel * taper
    out = np.clip(out, 0, None)
    s = out.sum()
    return out / s if s > 0 else kernel


def _kernel_shape(kernel):
    """FWHM_x/y, theta, elongation of a kernel from its second moments."""
    ny, nx = kernel.shape
    yy, xx = np.mgrid[:ny, :nx]
    w = np.clip(kernel, 0, None)
    tot = w.sum()
    cx, cy = (xx * w).sum() / tot, (yy * w).sum() / tot
    dx, dy = xx - cx, yy - cy
    ixx = (dx ** 2 * w).sum() / tot
    iyy = (dy ** 2 * w).sum() / tot
    ixy = (dx * dy * w).sum() / tot
    t, det = ixx + iyy, ixx * iyy - ixy ** 2
    disc = np.sqrt(max(0.25 * t ** 2 - det, 0))
    lam1, lam2 = 0.5 * t + disc, max(0.5 * t - disc, 1e-6)
    theta = 0.5 * np.arctan2(2 * ixy, ixx - iyy)
    f1, f2 = 2.355 * np.sqrt(lam1), 2.355 * np.sqrt(lam2)
    return f1, f2, float(theta), float(f1 / f2)


def _odd(n, lo=25):
    n = int(max(n, lo))
    return n if n % 2 == 1 else n + 1


# ─────────────────────────────────────────────────────────────────────────────
# Main API
# ─────────────────────────────────────────────────────────────────────────────

def measure_psf(
    image_data,
    saturation=None,
    seeing_hint_px=None,
    threshold_sigma=5.0,
    snr_min=20.0,
    max_stars=25,
    kernel_size=None,
    catalog_xy=None,
    catalog_match_rad=3.0,
    bkg_box=64,
):
    """
    Measure the PSF of an image.  Never raises for data reasons; always
    returns a normalised kernel.

    Parameters
    ----------
    image_data : 2-D ndarray
    saturation : float or None — ADU ceiling; default min(60000, 0.9*max).
    seeing_hint_px : float or None — header seeing in PIXELS (caller
        disambiguates units, e.g. via _fwhm_px_guess); used only as the
        last-resort fallback FWHM and as a consistency diagnostic.
    threshold_sigma : detection threshold in background-RMS units.
    snr_min : minimum peak/rms for PSF stars (PSFEx SAMPLE_MINSN ~ 20).
    max_stars : cap on stars used for the empirical stack.
    kernel_size : odd int or None — output kernel side; default 8*FWHM, >=25.
    catalog_xy : (N,2) array or None — pixel positions of KNOWN stars
        (e.g. PS1 matches).  When given, PSF stars must lie within
        catalog_match_rad px of a catalog star (kills galaxy contamination).
    bkg_box : background estimation box size.

    Returns
    -------
    dict:
      kernel      2-D ndarray, odd, sum = 1
      psf         alias of kernel (back-compat with psf_measure)
      fwhm, fwhm_x, fwhm_y, theta, elongation, beta
      mode        'hybrid' | 'moffat' | 'analytic' | 'fallback'
      n_stars     stars used in the final stack (0 for analytic/fallback)
      stars_x, stars_y
      diagnostics dict (n_detect, fwhm_frame, fwhm_scatter, resid_rms,
                        centroid_shift_max, header_seeing_consistent,
                        saturation_used, warnings [...])
    """
    warnings_list = []
    data = np.asarray(image_data, dtype=float)
    if data.ndim != 2 or min(data.shape) < 32:
        return _analytic_result(seeing_hint_px, kernel_size,
                                ['image too small or not 2-D'], mode='fallback')
    bad = ~np.isfinite(data)
    if bad.any():
        data = data.copy()
        data[bad] = np.nanmedian(data)

    # NOTE: when saturation is None we rely solely on the unit-free
    # flat-top test — inventing an ADU ceiling breaks on processed/scaled
    # images (e.g. NOT frames with star peaks at 7e5 in float units)

    # ── C.0 background, detection, frame FWHM bootstrap ─────────────────────
    bkg2d, rms2d = _estimate_background(data, box_size=bkg_box)
    data_sub = data - bkg2d
    pos_rms = rms2d[rms2d > 0]
    rms_scale = float(np.median(pos_rms)) if pos_rms.size else 1e-3
    # floor the rms so zero-padded alignment borders can't flood the
    # detection with noise peaks
    snr_img = data_sub / np.maximum(rms2d, 0.3 * rms_scale)
    # renormalise: on resampled images the noise is CORRELATED, so the
    # per-pixel rms underestimates it and every "SNR" is inflated — force
    # the background of the SNR image to unit robust scatter
    _bgpix = snr_img[np.isfinite(snr_img) & (np.abs(snr_img) < 10)]
    if _bgpix.size > 1000:
        _s0 = float(1.4826 * np.median(np.abs(_bgpix - np.median(_bgpix))))
        if _s0 > 1.2:
            snr_img = snr_img / _s0
            warnings_list.append(
                f'correlated noise: SNR scale deflated by {_s0:.2f}x')

    x, y, det_snr = _detect_sources(snr_img, threshold_sigma)
    n_detect = len(x)
    if n_detect == 0:
        warnings_list.append('no sources detected')
        return _analytic_result(seeing_hint_px, kernel_size, warnings_list,
                                mode='fallback')

    # frame FWHM from half-max areas of the brightest CR-free sources
    # (immune to cosmic rays and to Moffat-wing moment inflation)
    nbr = np.array([_neighbor_ratio(data_sub, x[i], y[i])
                    for i in range(n_detect)])
    boot_idx = np.where((det_snr > snr_min) & (nbr > 0.35))[0][:300]
    if len(boot_idx) < 5:
        boot_idx = np.where(nbr > 0.35)[0][:300]
    if len(boot_idx) < 5:
        boot_idx = np.arange(min(n_detect, 300))
    fwhm_hm = np.full(n_detect, np.nan)
    for i in boot_idx:
        fwhm_hm[i] = _fwhm_halfmax(data_sub, x[i], y[i])
    good_fw = fwhm_hm[np.isfinite(fwhm_hm) & (fwhm_hm > 1.3) & (fwhm_hm < 20.0)]
    if len(good_fw) == 0:
        warnings_list.append('no finite FWHM measurements')
        return _analytic_result(seeing_hint_px, kernel_size, warnings_list,
                                mode='fallback')
    # stars are the SMALLEST dense mode that VERIFIES as a consistent PSF —
    # in deep frames galaxies outnumber stars (median lands on galaxies),
    # while CR tracks can form an even lower mode (blind smallest fails).
    # Walk modes in ascending order; accept the first whose brightest
    # members Moffat-fit consistently above the physical 1.5 px floor.
    _ny0, _nx0 = data_sub.shape
    _modes = _mode_centers(good_fw, bin_w=0.4)
    fwhm_frame, fwhm_scatter = np.nan, np.nan
    for _mc in _modes:
        _win = max(0.8, 0.25 * _mc)
        _mem = [i for i in boot_idx
                if np.isfinite(fwhm_hm[i]) and abs(fwhm_hm[i] - _mc) <= _win]
        if len(_mem) < 3 and len(_modes) > 1:
            continue
        _fit_fw = []
        for i in _mem[:6]:
            _r = int(max(8, np.ceil(3.0 * _mc)))
            _xi, _yi = int(round(x[i])), int(round(y[i]))
            if not (_r <= _xi < _nx0 - _r and _r <= _yi < _ny0 - _r):
                continue
            _st = data_sub[_yi - _r: _yi + _r + 1, _xi - _r: _xi + _r + 1]
            _f = _fit_moffat(_st, max(_mc, 1.5))
            if (_f.get('ok') and _f['fwhm'] >= 1.5 and
                    np.hypot(_f['x0'] - _r, _f['y0'] - _r) <= 2.0):
                _fit_fw.append(_f['fwhm'])
        if len(_fit_fw) >= min(3, max(1, len(_mem))):
            _med = float(np.median(_fit_fw))
            _mad = 1.4826 * float(np.median(np.abs(np.array(_fit_fw) - _med)))
            if _med >= 1.5 and (len(_fit_fw) < 3 or _mad <= 0.35 * _med):
                fwhm_frame, fwhm_scatter = _med, _mad
                break
    if not np.isfinite(fwhm_frame):
        warnings_list.append('no mode verified as stellar — using smallest mode')
        fwhm_frame, fwhm_scatter = _smallest_mode(good_fw, bin_w=0.4)

    # candidate set: SNR window AND half-max FWHM within a band around the
    # frame value (rejects cosmic rays / galaxies before the size locus)
    cand0 = [i for i in range(n_detect)
             if det_snr[i] > snr_min and nbr[i] > 0.35
             and np.isfinite(fwhm_hm[i])
             and 0.55 * fwhm_frame <= fwhm_hm[i] <= 1.8 * fwhm_frame]
    cand0 = cand0[:300]

    # raw-peak (for saturation) + annulus-subtracted moments per candidate
    rad_mom = max(6.0, 3.0 * fwhm_frame)
    moms = []
    for i in cand0:
        m = _source_moments(data_sub, x[i], y[i], rad_mom)
        if m is None:
            continue
        xi, yi = int(round(x[i])), int(round(y[i]))
        raw_peak = float(np.nanmax(data[max(0, yi - 1): yi + 2,
                                        max(0, xi - 1): xi + 2]))
        m.update(x=x[i], y=y[i], peak=raw_peak, snr=det_snr[i],
                 fwhm_hm=fwhm_hm[i])
        moms.append(m)
    if len(moms) == 0:
        warnings_list.append('no measurable candidates')
        return _analytic_result(fwhm_frame, kernel_size, warnings_list,
                                mode='analytic')
    snr_all = np.array([m['snr'] for m in moms])

    header_consistent = None
    if seeing_hint_px is not None and np.isfinite(seeing_hint_px):
        header_consistent = bool(0.5 <= seeing_hint_px / max(fwhm_frame, 1e-3) <= 2.0)
        if not header_consistent:
            warnings_list.append(
                f'header seeing {seeing_hint_px:.1f}px inconsistent with '
                f'measured {fwhm_frame:.1f}px')

    # ── C.1 star selection ──────────────────────────────────────────────────
    stamp_size = _odd(8 * fwhm_frame, lo=21)
    half = stamp_size // 2
    ny, nx = data.shape

    rh_arr = np.array([m.get('rh', np.nan) for m in moms])
    mag_ok = np.isfinite(rh_arr) & (snr_all > snr_min)
    cand_idx = np.where(mag_ok)[0]
    # size locus: unresolved sources form the SMALLEST half-light-radius
    # mode (PSFEx/LSST objectSize); galaxies sit above it, so use the low
    # mode, never a median over a galaxy-rich candidate list
    rh_locus, rh_mad = _smallest_mode(rh_arr[cand_idx], bin_w=0.25)
    tol = (max(0.2 * rh_locus, 3.0 * rh_mad, 0.35)
           if np.isfinite(rh_locus) else np.inf)

    # elongation statistics from LOCUS MEMBERS only (stars) — the full
    # candidate list is galaxy-contaminated on deep frames
    _locus_el = np.array([moms[i]['elongation'] for i in cand_idx
                          if np.isfinite(rh_locus)
                          and abs(moms[i]['rh'] - rh_locus) <= tol
                          and np.isfinite(moms[i]['elongation'])])
    if len(_locus_el) >= 3:
        el_med = float(np.median(_locus_el))
        el_mad = 1.4826 * float(np.median(np.abs(_locus_el - el_med)))
    else:
        el_med, el_mad = 1.2, 0.3
    _locus_th = np.array([moms[i]['theta'] for i in cand_idx
                          if np.isfinite(rh_locus)
                          and abs(moms[i]['rh'] - rh_locus) <= tol
                          and np.isfinite(moms[i]['theta'])])
    frame_elong = el_med
    frame_theta = float(np.median(_locus_th)) if len(_locus_th) else 0.0

    base = []
    for i in cand_idx:
        m = moms[i]
        if np.isfinite(rh_locus) and abs(m['rh'] - rh_locus) > tol:
            continue                              # off the stellar locus
        if saturation is not None and m['peak'] > 0.9 * saturation:
            continue                              # nonlinearity margin
        # flat-top saturation test needs noise << 2% of peak to discriminate;
        # below SNR 100 it fires on noise ties (and such stars can't saturate)
        if m['snr'] > 100 and _is_flattop(data, m['x'], m['y'], fwhm_frame):
            continue                              # saturated plateau (unit-free)
        if not (half + 2 <= m['x'] < nx - half - 2 and
                half + 2 <= m['y'] < ny - half - 2):
            continue                              # edge margin
        if np.isfinite(m['elongation']) and \
                abs(m['elongation'] - el_med) > max(3.0 * el_mad, 0.4):
            continue                              # shape outlier vs frame
        if catalog_xy is not None and len(catalog_xy) > 0:
            cd = np.hypot(catalog_xy[:, 0] - m['x'], catalog_xy[:, 1] - m['y'])
            if cd.min() > catalog_match_rad:
                continue                          # not a known catalog star
        base.append(m)

    # graduated isolation: start strict; relax on crowded/deep frames until
    # enough stars survive (off-centre neighbours in the outer stamp are
    # suppressed by the pixelwise sigma-clipped stack anyway)
    selected = []
    for rad_f, frac in ((4.0, 0.10), (2.5, 0.20), (1.5, 0.50)):
        selected = []
        for m in base:
            d = np.hypot(x - m['x'], y - m['y'])
            near = (d > 0.5) & (d < rad_f * fwhm_frame)
            if np.any(det_snr[near] > np.maximum(8.0, frac * m['snr'])):
                continue
            selected.append(m)
        if len(selected) >= 5:
            break
    if selected and (rad_f, frac) != (4.0, 0.10):
        warnings_list.append(
            f'crowded field: isolation relaxed to {rad_f:.1f}xFWHM')
    selected.sort(key=lambda m: m['peak'], reverse=True)
    selected = selected[:max_stars]

    if len(selected) < 2:
        warnings_list.append(f'only {len(selected)} PSF stars — analytic mode')
        return _analytic_result(
            fwhm_frame, kernel_size, warnings_list, mode='analytic',
            elongation=frame_elong, theta=frame_theta,
            extra=dict(n_detect=n_detect, fwhm_frame=fwhm_frame,
                       fwhm_scatter=fwhm_scatter,
                       header_seeing_consistent=header_consistent,
                       saturation_used=saturation))

    # ── C.2 stage 1: Moffat-centroided sigma-clipped stack ──────────────────
    stamps, fits_ok, cshift = [], [], []
    for m in selected:
        xi, yi = int(round(m['x'])), int(round(m['y']))
        st = data_sub[yi - half: yi + half + 1, xi - half: xi + half + 1].astype(float)
        if st.shape != (stamp_size, stamp_size):
            continue
        yy, xx = np.mgrid[:stamp_size, :stamp_size]
        rr = np.hypot(xx - half, yy - half)
        ann = st[(rr > 3.0 * fwhm_frame) & (rr <= min(4.0 * fwhm_frame, half))]
        if ann.size >= 8:
            st = st - float(np.median(ann))
        fit = _fit_moffat(st, fwhm_frame)
        if not fit.get('ok'):
            continue
        dxc, dyc = fit['x0'] - half, fit['y0'] - half
        shift_mag = float(np.hypot(dxc, dyc))
        if shift_mag > 1.5:
            continue                              # centroid ran away — drop
        if not (0.5 * fwhm_frame <= fit['fwhm'] <= 2.0 * fwhm_frame):
            continue                              # not the frame's PSF (CR/blend)
        shifted = ndimage.shift(st, (-dyc, -dxc), order=3, mode='constant', cval=0.0)
        norm = np.clip(shifted, 0, None)[rr <= 3.0 * fwhm_frame].sum()
        if norm <= 0:
            continue
        stamps.append(shifted / norm)
        fits_ok.append(fit)
        cshift.append(shift_mag)

    if len(stamps) < 2:
        warnings_list.append('star fits failed — analytic mode')
        return _analytic_result(
            fwhm_frame, kernel_size, warnings_list, mode='analytic',
            elongation=frame_elong, theta=frame_theta,
            extra=dict(n_detect=n_detect, fwhm_frame=fwhm_frame,
                       fwhm_scatter=fwhm_scatter,
                       header_seeing_consistent=header_consistent,
                       saturation_used=saturation))

    stamps = np.array(stamps)
    # Piff-style loop: stack, reject worst star by chi2, restack (<=3 rounds)
    for _ in range(3):
        clipped = sigma_clip(stamps, sigma=3.0, axis=0, maxiters=3)
        stack = np.ma.median(clipped, axis=0).filled(0.0)
        if len(stamps) <= 3:
            break
        chi = np.array([np.mean((s - stack) ** 2) for s in stamps])
        worst = int(np.argmax(chi))
        if chi[worst] > 4.0 * np.median(chi):
            stamps = np.delete(stamps, worst, axis=0)
            fits_ok.pop(worst)
            cshift.pop(worst)
        else:
            break
    n_used = len(stamps)
    resid_rms = float(np.median([np.sqrt(np.mean((s - stack) ** 2)) for s in stamps])
                      / max(stack.max(), 1e-9))
    per_star_fwhm = np.array([f['fwhm'] for f in fits_ok])
    star_fwhm_scatter = (float(1.4826 * np.median(
        np.abs(per_star_fwhm - np.median(per_star_fwhm))))
        if n_used > 2 else 0.0)

    # ── C.3 stage 2: Moffat fit to the stack ────────────────────────────────
    per_star_med = float(np.median([f['fwhm'] for f in fits_ok]))
    stack_fit = _fit_moffat(stack, fwhm_frame)
    if stack_fit.get('ok'):
        geo = float(np.sqrt(stack_fit['fwhm_x'] * stack_fit['fwhm_y']))
        # sanity: the stack fit must agree with the per-star consensus —
        # a collapsed/expanded fit means the stack was contaminated
        if not (1.2 <= geo <= 25.0) or \
                abs(geo - per_star_med) > 0.5 * per_star_med:
            stack_fit = dict(ok=False)
            warnings_list.append(
                f'stack fit FWHM {geo:.1f}px inconsistent with per-star '
                f'median {per_star_med:.1f}px — using per-star medians')
    if not stack_fit.get('ok'):
        stack_fit = dict(
            ok=True, beta=float(np.median([f['beta'] for f in fits_ok])),
            fwhm_x=float(np.median([f['fwhm_x'] for f in fits_ok])),
            fwhm_y=float(np.median([f['fwhm_y'] for f in fits_ok])),
            theta=float(np.median([f['theta'] for f in fits_ok])),
            x0=float(half), y0=float(half), rms=np.nan)
        if 'stack fit' not in ' '.join(warnings_list):
            warnings_list.append('stack Moffat fit failed; using per-star medians')
        geo = float(np.sqrt(stack_fit['fwhm_x'] * stack_fit['fwhm_y']))
        if not (1.2 <= geo <= 25.0):
            warnings_list.append(
                f'per-star FWHM {geo:.1f}px unphysical — analytic mode')
            return _analytic_result(
                fwhm_frame, kernel_size, warnings_list, mode='analytic',
                elongation=frame_elong, theta=frame_theta,
                extra=dict(n_detect=n_detect, fwhm_frame=fwhm_frame,
                           fwhm_scatter=fwhm_scatter,
                           header_seeing_consistent=header_consistent,
                           saturation_used=saturation))

    ksize = _odd(kernel_size if kernel_size else 8 * fwhm_frame)
    moffat_k = _moffat_kernel(ksize, stack_fit['fwhm_x'], stack_fit['fwhm_y'],
                              stack_fit['theta'], stack_fit['beta'])

    # ── C.4 stage 3: hybrid = analytic + capped residual ────────────────────
    mode = 'moffat'
    kernel = moffat_k
    if n_used >= 5 and np.isfinite(resid_rms) and resid_rms < 0.25:
        stack_n = stack.copy()
        s = np.clip(stack_n, 0, None).sum()
        stack_n = stack_n / s if s > 0 else stack_n
        # put the stack onto the kernel grid (both odd, centre-aligned)
        stack_k = np.zeros_like(moffat_k)
        hs, hk = stamp_size // 2, ksize // 2
        h = min(hs, hk)
        stack_k[hk - h: hk + h + 1, hk - h: hk + h + 1] = \
            stack_n[hs - h: hs + h + 1, hs - h: hs + h + 1]
        resid = stack_k - moffat_k
        resid = ndimage.gaussian_filter(resid, sigma=0.7)
        cap = 0.3 * moffat_k.max()
        resid = np.clip(resid, -cap, cap)
        kernel = moffat_k + resid
        mode = 'hybrid'
    elif n_used >= 5:
        warnings_list.append(
            f'residual rms {resid_rms:.2f} too high — Moffat-only kernel')

    fwhm_geo = float(np.sqrt(stack_fit['fwhm_x'] * stack_fit['fwhm_y']))
    kernel = _apodize(kernel, fwhm_geo)
    kf1, kf2, ktheta, kelong = _kernel_shape(kernel)

    diagnostics = dict(
        mode=mode, n_detect=n_detect, n_stars_used=n_used,
        fwhm_frame=fwhm_frame, fwhm_scatter=fwhm_scatter,
        star_fwhm_scatter=star_fwhm_scatter, resid_rms=resid_rms,
        centroid_shift_max=float(max(cshift)) if cshift else np.nan,
        moffat_beta=float(stack_fit['beta']),
        header_seeing_consistent=header_consistent,
        saturation_used=float(saturation) if saturation is not None else None,
        catalog_matched=catalog_xy is not None,
        warnings=warnings_list,
    )
    return dict(
        kernel=kernel, psf=kernel,
        fwhm=fwhm_geo,
        fwhm_x=float(stack_fit['fwhm_x']), fwhm_y=float(stack_fit['fwhm_y']),
        theta=float(stack_fit['theta']),
        elongation=float(max(stack_fit['fwhm_x'], stack_fit['fwhm_y']) /
                         max(min(stack_fit['fwhm_x'], stack_fit['fwhm_y']), 1e-6)),
        beta=float(stack_fit['beta']),
        mode=mode, n_stars=n_used,
        stars_x=np.array([f['x0'] for f in fits_ok]),
        stars_y=np.array([f['y0'] for f in fits_ok]),
        diagnostics=diagnostics,
    )


def _analytic_result(fwhm_px, kernel_size, warnings_list, mode='analytic',
                     elongation=1.0, theta=0.0, extra=None):
    """Build the analytic/fallback return: elliptical (or circular) Moffat
    from whatever shape information is available."""
    if fwhm_px is None or not np.isfinite(fwhm_px):
        fwhm_px = 4.0
        warnings_list = warnings_list + ['no seeing information; FWHM=4px assumed']
    fwhm_px = float(np.clip(fwhm_px, FWHM_MIN_PX, FWHM_MAX_PX))
    elongation = float(np.clip(elongation, 1.0, 3.0)) if np.isfinite(elongation) else 1.0
    # split geometric-mean FWHM into major/minor at the given elongation
    fwhm_maj = fwhm_px * np.sqrt(elongation)
    fwhm_min = fwhm_px / np.sqrt(elongation)
    ksize = _odd(kernel_size if kernel_size else 8 * fwhm_px)
    kernel = _moffat_kernel(ksize, fwhm_maj, fwhm_min,
                            theta if np.isfinite(theta) else 0.0,
                            _MOFFAT_BETA_DEFAULT)
    diagnostics = dict(mode=mode, n_stars_used=0, warnings=warnings_list,
                       moffat_beta=_MOFFAT_BETA_DEFAULT)
    if extra:
        diagnostics.update(extra)
    return dict(kernel=kernel, psf=kernel, fwhm=fwhm_px,
                fwhm_x=float(fwhm_maj), fwhm_y=float(fwhm_min),
                theta=float(theta if np.isfinite(theta) else 0.0),
                elongation=elongation, beta=_MOFFAT_BETA_DEFAULT,
                mode=mode, n_stars=0,
                stars_x=np.array([]), stars_y=np.array([]),
                diagnostics=diagnostics)


def measure_psf_from_fits(filepath, ext=0, **kwargs):
    """Load a FITS image and measure its PSF (see measure_psf)."""
    with fits.open(filepath) as hdul:
        data = hdul[ext].data
        header = hdul[ext].header
    if data is None:
        for hdu in hdul[1:]:
            if getattr(hdu, 'data', None) is not None:
                data = hdu.data
                header = hdu.header
                break
    if data is None:
        raise ValueError(f'no image data in {filepath}')
    if 'saturation' not in kwargs:
        sat = header.get('SATURATE', header.get('SATURLEV'))
        if sat:
            kwargs['saturation'] = 0.95 * float(sat)
    return measure_psf(np.asarray(data, float), **kwargs)


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        res = measure_psf_from_fits(sys.argv[1])
    else:
        rng = np.random.default_rng(7)
        ny = nx = 600
        img = rng.poisson(150.0, (ny, nx)).astype(float)
        yy, xx = np.mgrid[:ny, :nx]
        beta_t, fw_t = 3.0, 4.6
        a_t = _alpha_from_fwhm(fw_t, beta_t)
        for _ in range(25):
            x0, y0 = rng.integers(50, nx - 50, 2)
            amp = rng.uniform(2000, 30000)
            img += amp * (1 + ((xx - x0) ** 2 + ((yy - y0) * 1.15) ** 2) / a_t ** 2) ** -beta_t
        res = measure_psf(img, saturation=60000)
        print(f'true FWHM ~{fw_t:.2f}px (y squashed 1.15x)')
    d = res['diagnostics']
    print(f"mode={res['mode']}  FWHM={res['fwhm']:.2f}px "
          f"({res['fwhm_x']:.2f} x {res['fwhm_y']:.2f})  "
          f"elong={res['elongation']:.2f}  beta={res['beta']:.2f}  "
          f"n_stars={res['n_stars']}")
    print('diagnostics:', {k: v for k, v in d.items() if k != 'warnings'})
    if d['warnings']:
        print('warnings:', d['warnings'])
