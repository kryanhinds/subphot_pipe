#!/usr/bin/env python3
"""
Astrometric recalibration using detected image sources and Pan-STARRS1 catalog stars.

Workflow:
1) Detect stars in image (photutils DAOStarFinder).
2) Query PS1 around image center (astroquery.mast Catalogs).
3) Coarse match in sky coordinates using initial WCS, estimate bulk offset.
4) Re-match with tighter tolerance.
5) Fit WCS from matched (x, y) <-> (RA, Dec) points using astropy.
"""

import argparse
import os
import sys
import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.io.votable import parse_single_table
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import fit_wcs_from_points

from scipy.ndimage import maximum_filter, label, center_of_mass
import requests

try:
    from astroquery.mast import Catalogs
    HAS_ASTROQUERY = True
except Exception:
    HAS_ASTROQUERY = False

try:
    from photutils.detection import DAOStarFinder
    HAS_PHOTUTILS = True
except Exception:
    HAS_PHOTUTILS = False


def _float_header_value(header, keys, default=None):
    for key in keys:
        if key in header:
            try:
                return float(header[key])
            except Exception:
                pass
    return default


def _safe_float(value):
    if value is None:
        return None
    s = str(value).strip()
    if s == "" or s.lower() in {"none", "null", "nan"}:
        return None
    try:
        return float(s)
    except Exception:
        return None


def _estimate_seeing_pix(header, pixscale_arcsec):
    pix = _float_header_value(header, ["SEEINGPIX", "FWHMPIX", "L1SEEPX"], None)
    if pix and pix > 0:
        return pix

    sec = _float_header_value(
        header,
        ["L1SEESEC", "SEEING", "FWHM", "ESO TEL AMBI FWHM END"],
        None,
    )
    if sec and sec > 0 and pixscale_arcsec > 0:
        return sec / pixscale_arcsec

    return 2.8


def _find_colname(table, preferred):
    cols = {c.lower(): c for c in table.colnames}
    for name in preferred:
        if name.lower() in cols:
            return cols[name.lower()]
    return None


def _column_to_float_array(col):
    """
    Convert an astropy column (possibly masked/int/string) to float ndarray with NaN for invalid entries.
    """
    out = np.full(len(col), np.nan, dtype=float)
    if np.ma.isMaskedArray(col):
        data = np.asarray(col.data)
        mask = np.ma.getmaskarray(col)
    else:
        data = np.asarray(col)
        mask = np.zeros(len(data), dtype=bool)

    for i, val in enumerate(data):
        if mask[i]:
            continue
        fval = _safe_float(val)
        if fval is not None:
            out[i] = fval
    return out


def panstarrs_query(
    ra_deg,
    dec_deg,
    rad_deg,
    logger=None,
    mindet=1,
    maxsources=10000,
    server="https://archive.stsci.edu/panstarrs/search.php",
    cache_dir="ps_catalogs",
):
    """
    Query Pan-STARRS @ MAST via VOTable search endpoint and cache XML locally.
    Returns astropy.table.Table.
    """
    os.makedirs(cache_dir, exist_ok=True)
    cache_name = f"ps_{ra_deg:.6f}_{dec_deg:.6f}_{rad_deg:.5f}.xml"
    cache_path = os.path.join(cache_dir, cache_name)

    if not os.path.exists(cache_path):
        r = requests.get(
            server,
            params={
                "RA": ra_deg,
                "DEC": dec_deg,
                "SR": rad_deg,
                "max_records": maxsources,
                "outputformat": "VOTable",
                "ndetections": f">{mindet}",
            },
            timeout=90,
        )
        r.raise_for_status()
        with open(cache_path, "w", encoding="utf-8") as outf:
            outf.write(r.text)
        if logger is not None:
            logger.info("PS1 catalog downloaded: %s", cache_path)
    else:
        if logger is not None:
            logger.info("PS1 catalog cached: %s", cache_path)

    data = parse_single_table(cache_path)
    return data.to_table(use_names_over_ids=True)


def detect_sources(data, header, saturation, threshold_sigma=5.0, max_sources=600):
    if data is None:
        raise RuntimeError("No image data found in FITS.")

    mean, median, std = sigma_clipped_stats(data, sigma=3.0, maxiters=5)
    pixscale = _float_header_value(header, ["PIXSCALE"], 0.394)
    seeing_pix = _estimate_seeing_pix(header, pixscale)
    fwhm = float(np.clip(seeing_pix, 1.5, 8.0))

    def detect_stars_with_size(image, threshold_sigma=5.0, fwhm=3.0, border_mask=None, n_max=None, avoid_radius=None):
        """
        Detect stars and estimate approximate size while suppressing overlaps.
        """
        finite = np.isfinite(image)
        if not np.any(finite):
            return np.array([]), np.array([]), np.array([])
        mean = np.nanmean(image)
        median = np.nanmedian(image)
        std = np.nanstd(image)
        if not np.isfinite(std) or std <= 0:
            return np.array([]), np.array([]), np.array([])
        threshold = median + threshold_sigma * std

        img = np.array(image, dtype=float, copy=True)
        img[~finite] = median

        if border_mask is None:
            valid_mask = np.ones_like(image, dtype=bool)
        else:
            valid_mask = border_mask.copy()
        valid_mask &= finite

        xy_list = []
        flux_list = []
        size_list = []

        while True:
            size_filter = max(2, int(np.ceil(fwhm)))
            work = img * valid_mask
            local_max = maximum_filter(work, size=size_filter) == work
            detect = local_max & (img > threshold) & valid_mask

            if saturation and saturation > 0:
                detect &= img < saturation

            labeled, n_labels = label(detect)
            if n_labels == 0:
                break

            peak_flux = 0.0
            peak_idx = -1
            for region_id in range(1, n_labels + 1):
                mask = labeled == region_id
                total_flux = float(np.sum(img[mask]))
                if total_flux > peak_flux:
                    peak_flux = total_flux
                    peak_idx = region_id

            if peak_idx == -1:
                break

            mask = labeled == peak_idx
            try:
                y_c, x_c = center_of_mass(img, labels=mask, index=1)
            except Exception:
                ys, xs = np.nonzero(mask)
                if len(xs) == 0:
                    valid_mask[mask] = False
                    continue
                x_c = float(np.mean(xs))
                y_c = float(np.mean(ys))

            ys, xs = np.nonzero(mask)
            fluxes = img[ys, xs]
            flux_sum = np.sum(fluxes)
            if not np.isfinite(flux_sum) or flux_sum <= 0:
                valid_mask[mask] = False
                continue

            x_mean = float(np.sum(xs * fluxes) / flux_sum)
            y_mean = float(np.sum(ys * fluxes) / flux_sum)
            sigma_x = float(np.sqrt(np.sum(fluxes * (xs - x_mean) ** 2) / flux_sum))
            sigma_y = float(np.sqrt(np.sum(fluxes * (ys - y_mean) ** 2) / flux_sum))
            star_size = float(np.mean([sigma_x, sigma_y]))

            xy_list.append([x_c, y_c])
            flux_list.append(float(np.sum(img[mask])))
            size_list.append(star_size)

            if avoid_radius is not None:
                y_grid, x_grid = np.indices(image.shape)
                dist2 = (x_grid - x_c) ** 2 + (y_grid - y_c) ** 2
                mask_radius = dist2 <= (avoid_radius + star_size) ** 2
                valid_mask[mask_radius] = False
            else:
                valid_mask[mask] = False

            if n_max is not None and len(xy_list) >= n_max:
                break

        return np.array(xy_list), np.array(flux_list), np.array(size_list)

    ny, nx = data.shape
    # SEDM frames often have unstable edge structure; keep detections well inside frame.
    border = max(40, int(5 * fwhm))
    border_mask = np.zeros_like(data, dtype=bool)
    border_mask[border:ny - border, border:nx - border] = True

    xy, flux, size = detect_stars_with_size(
        data,
        threshold_sigma=threshold_sigma,
        fwhm=fwhm,
        border_mask=border_mask,
        n_max=max_sources,
        avoid_radius=max(2.0, 0.8 * fwhm),
    )
    if len(xy) > 0:
        good = np.isfinite(flux) & np.isfinite(size)
        xy = xy[good]
        flux = flux[good]
        if len(xy) > 0:
            order = np.argsort(flux)[::-1][:max_sources]
            return xy[order, 0], xy[order, 1]

    # Fallback to DAOStarFinder if custom detector finds nothing.
    if HAS_PHOTUTILS:
        daofind = DAOStarFinder(
            fwhm=fwhm,
            threshold=threshold_sigma * std,
            exclude_border=True,
        )
        src = daofind(data - median)
        if src is not None and len(src) > 0:
            tab = Table(src)
            if "peak" in tab.colnames and saturation and saturation > 0:
                tab = tab[tab["peak"] < saturation]
            if len(tab) > 0:
                tab.sort("flux")
                tab.reverse()
                tab = tab[:max_sources]
                return np.array(tab["xcentroid"]), np.array(tab["ycentroid"])

    # Last-resort simple local-max detector with adaptive thresholds.
    finite = np.isfinite(data)
    if not np.any(finite):
        return np.array([]), np.array([])
    med = np.nanmedian(data)
    sig = np.nanstd(data)
    img = np.array(data, dtype=float, copy=True)
    img[~finite] = med
    yy, xx = np.indices(img.shape)
    valid = finite & border_mask
    for tscale in (threshold_sigma, 4.0, 3.0, 2.5):
        thr = med + tscale * sig
        peaks = (img == maximum_filter(img, size=max(3, int(np.ceil(fwhm))))) & (img > thr) & valid
        if saturation and saturation > 0:
            peaks &= img < saturation
        y, x = np.where(peaks)
        if len(x) == 0:
            continue
        f = img[y, x]
        order = np.argsort(f)[::-1][:max_sources]
        return x[order].astype(float), y[order].astype(float)

    return np.array([]), np.array([])


def query_ps1(center, radius_deg, band, mag_limit=21.5):
    mag_col = {
        "g": "gMeanPSFMag",
        "r": "rMeanPSFMag",
        "i": "iMeanPSFMag",
        "z": "zMeanPSFMag",
        "y": "yMeanPSFMag",
    }.get(band.lower(), "rMeanPSFMag")

    # Primary path: PS1 VOTable endpoint with local caching.
    try:
        tab = panstarrs_query(
            float(center.ra.deg),
            float(center.dec.deg),
            float(radius_deg),
            mindet=2,
            maxsources=10000,
        )
    except Exception:
        tab = None

    # Fallback path: astroquery if available.
    if (tab is None or len(tab) == 0) and HAS_ASTROQUERY:
        tab = Catalogs.query_region(
            center,
            radius=radius_deg * u.deg,
            catalog="Panstarrs",
            data_release="dr2",
            table="mean",
        )

    if tab is None or len(tab) == 0:
        return SkyCoord([], [], unit="deg")

    ra_col = _find_colname(tab, ["raMean", "ra"])
    dec_col = _find_colname(tab, ["decMean", "dec"])
    m_col = _find_colname(tab, [mag_col, "rMeanPSFMag", "rMeanKronMag", "rmag"])
    n_col = _find_colname(tab, ["nDetections"])

    if ra_col is None or dec_col is None:
        return SkyCoord([], [], unit="deg")

    ra_arr = _column_to_float_array(tab[ra_col])
    dec_arr = _column_to_float_array(tab[dec_col])
    keep = np.isfinite(ra_arr) & np.isfinite(dec_arr)

    if m_col is not None:
        m_arr = _column_to_float_array(tab[m_col])
        keep &= np.isfinite(m_arr)
        keep &= m_arr < mag_limit
        keep &= m_arr > 12.0

    if n_col is not None:
        n_arr = _column_to_float_array(tab[n_col])
        keep &= np.isfinite(n_arr)
        keep &= n_arr > 1

    if np.count_nonzero(keep) == 0:
        return SkyCoord([], [], unit="deg")

    return SkyCoord(ra_arr[keep], dec_arr[keep], unit="deg", frame="icrs")


def coarse_offset(det_world, cat_world, coarse_arcsec=120.0):
    idx, sep2d, _ = det_world.match_to_catalog_sky(cat_world)
    keep = sep2d.arcsec < coarse_arcsec
    if np.count_nonzero(keep) < 6:
        return 0.0, 0.0, keep

    # Keep dRA in true degree units (no cos(dec) factor) so it can be directly applied.
    dra = (cat_world[idx[keep]].ra.deg - det_world[keep].ra.deg)
    ddec = cat_world[idx[keep]].dec.deg - det_world[keep].dec.deg

    med_dra = float(np.median(dra))
    med_ddec = float(np.median(ddec))
    return med_dra, med_ddec, keep


def _estimate_offset_histogram(det_world, cat_world, center, n_det=120, n_cat=300, bin_arcsec=3.0):
    """
    Estimate bulk sky offset in tangent-plane arcsec using a 2D difference histogram.
    """
    n_det = min(n_det, len(det_world))
    n_cat = min(n_cat, len(cat_world))
    if n_det < 5 or n_cat < 5:
        return 0.0, 0.0

    cdet = det_world[:n_det]
    ccat = cat_world[:n_cat]
    cosd0 = np.cos(np.deg2rad(center.dec.deg))

    dx_det = (cdet.ra.deg - center.ra.deg) * cosd0 * 3600.0
    dy_det = (cdet.dec.deg - center.dec.deg) * 3600.0
    dx_cat = (ccat.ra.deg - center.ra.deg) * cosd0 * 3600.0
    dy_cat = (ccat.dec.deg - center.dec.deg) * 3600.0

    ddx = (dx_cat[:, None] - dx_det[None, :]).ravel()
    ddy = (dy_cat[:, None] - dy_det[None, :]).ravel()

    lim = 900.0
    keep = np.isfinite(ddx) & np.isfinite(ddy) & (np.abs(ddx) < lim) & (np.abs(ddy) < lim)
    ddx = ddx[keep]
    ddy = ddy[keep]
    if len(ddx) == 0:
        return 0.0, 0.0

    nb = int((2 * lim) / bin_arcsec) + 1
    h, xed, yed = np.histogram2d(ddx, ddy, bins=[nb, nb], range=[[-lim, lim], [-lim, lim]])
    ip = np.unravel_index(np.argmax(h), h.shape)
    offx = 0.5 * (xed[ip[0]] + xed[ip[0] + 1])
    offy = 0.5 * (yed[ip[1]] + yed[ip[1] + 1])

    dra_deg = offx / (3600.0 * cosd0)
    ddec_deg = offy / 3600.0
    return float(dra_deg), float(ddec_deg)


def apply_offset(coords, dra_deg, ddec_deg):
    ra = (coords.ra.deg + dra_deg) % 360.0
    dec = coords.dec.deg + ddec_deg
    return SkyCoord(ra, dec, unit="deg", frame="icrs")


def solve_one(
    filename,
    output=None,
    band="r",
    mag_limit=21.5,
    threshold_sigma=5.0,
    coarse_arcsec=120.0,
    fine_arcsec=2.2,
    min_matches=12,
    saturation=None,
):
    hdul = fits.open(filename)
    hdu = hdul[0]
    data = hdu.data
    header = hdu.header.copy()

    if data is not None and data.ndim > 2:
        data = np.squeeze(data)
    if data is None or data.ndim != 2:
        hdul.close()
        raise RuntimeError(f"{filename}: expected 2D image data.")

    w0 = WCS(header)
    if not w0.has_celestial:
        hdul.close()
        raise RuntimeError(f"{filename}: no celestial WCS in header.")

    xc, yc = data.shape[1] / 2.0, data.shape[0] / 2.0
    c0 = w0.pixel_to_world(xc, yc)
    c0 = SkyCoord(c0.ra, c0.dec, frame="icrs")

    # Estimate query radius from image diagonal using current WCS pixel scale.
    pix_scales = w0.proj_plane_pixel_scales()
    ps_vals = []
    for p in pix_scales:
        try:
            ps_vals.append(float(np.abs(p.to(u.deg).value)))
        except Exception:
            try:
                ps_vals.append(float(np.abs(p)))
            except Exception:
                pass
    pscale = float(np.nanmedian(ps_vals) * 3600.0) if len(ps_vals) else np.nan
    if not np.isfinite(pscale) or pscale <= 0:
        pscale = 0.394
    diag_pix = np.hypot(data.shape[0], data.shape[1])
    radius_deg = max(0.08, (diag_pix * pscale / 3600.0) * 0.7)

    sat = saturation
    if sat is None:
        sat = _float_header_value(header, ["SATURATE", "SATURAT", "SATLEVEL", "SATUR_LEVEL"], 50000.0)

    x, y = detect_sources(data, header, sat, threshold_sigma=threshold_sigma)
    if len(x) < min_matches:
        hdul.close()
        raise RuntimeError(f"{filename}: only {len(x)} detections found.")
    det_world = w0.pixel_to_world(x, y)
    det_world = SkyCoord(det_world.ra, det_world.dec, frame="icrs")

    cat_world = query_ps1(c0, radius_deg, band=band, mag_limit=mag_limit)
    if len(cat_world) < min_matches:
        hdul.close()
        raise RuntimeError(f"{filename}: PS1 query returned too few stars ({len(cat_world)}).")

    dra, ddec, coarse_keep = coarse_offset(det_world, cat_world, coarse_arcsec=coarse_arcsec)
    shifted = apply_offset(det_world, dra, ddec)

    idx, sep2d, _ = shifted.match_to_catalog_sky(cat_world)
    good = sep2d.arcsec < fine_arcsec

    # Adaptive retries for difficult fields.
    used_fine_arcsec = fine_arcsec
    if np.count_nonzero(good) < min_matches:
        for trial in (max(fine_arcsec, 3.0), 4.0, 5.0):
            good_t = sep2d.arcsec < trial
            if np.count_nonzero(good_t) > np.count_nonzero(good):
                good = good_t
                used_fine_arcsec = trial
            if np.count_nonzero(good) >= min_matches:
                break

    # If still failing, estimate bulk offset by 2D histogram in tangent plane.
    if np.count_nonzero(good) < min_matches:
        dra2, ddec2 = _estimate_offset_histogram(det_world, cat_world, c0)
        shifted2 = apply_offset(det_world, dra2, ddec2)
        idx2, sep2d_2, _ = shifted2.match_to_catalog_sky(cat_world)
        best_good = good
        best_idx = idx
        best_sep = sep2d
        for trial in (fine_arcsec, 3.0, 4.0, 5.0):
            good_t = sep2d_2.arcsec < trial
            if np.count_nonzero(good_t) > np.count_nonzero(best_good):
                best_good = good_t
                best_idx = idx2
                best_sep = sep2d_2
                dra, ddec = dra2, ddec2
                used_fine_arcsec = trial
        idx, sep2d, good = best_idx, best_sep, best_good

    hard_min = min(8, min_matches)
    if np.count_nonzero(good) < hard_min:
        hdul.close()
        raise RuntimeError(
            f"{filename}: only {np.count_nonzero(good)} refined matches "
            f"(need >= {min_matches}, hard floor {hard_min})."
        )

    xg = x[good]
    yg = y[good]
    skyg = cat_world[idx[good]]

    try:
        w_new = fit_wcs_from_points(
            (xg, yg),
            skyg,
            projection=w0.celestial,
            sip_degree=3,
        )
    except Exception:
        w_new = fit_wcs_from_points((xg, yg), skyg, projection=w0.celestial)

    # Robust inlier refinement after first fit.
    for _ in range(2):
        pred = w_new.pixel_to_world(xg, yg)
        pred = SkyCoord(pred.ra, pred.dec, frame="icrs")
        sep = pred.separation(skyg).arcsec
        if len(sep) < hard_min:
            break
        med = np.median(sep)
        mad = np.median(np.abs(sep - med))
        clip = max(1.2, med + 3.0 * (1.4826 * mad))
        inl = sep < clip
        if np.count_nonzero(inl) < hard_min or np.count_nonzero(inl) == len(sep):
            break
        xg = xg[inl]
        yg = yg[inl]
        skyg = skyg[inl]
        try:
            w_new = fit_wcs_from_points(
                (xg, yg),
                skyg,
                projection=w0.celestial,
                sip_degree=3,
            )
        except Exception:
            w_new = fit_wcs_from_points((xg, yg), skyg, projection=w0.celestial)

    pred = w_new.pixel_to_world(xg, yg)
    pred = SkyCoord(pred.ra, pred.dec, frame="icrs")
    fit_sep = pred.separation(skyg).arcsec
    fit_rms = float(np.sqrt(np.mean(fit_sep**2))) if len(fit_sep) else 999.0
    if len(xg) < hard_min or fit_rms > 2.0:
        hdul.close()
        raise RuntimeError(
            f"{filename}: poor final fit quality (inliers={len(xg)}, rms={fit_rms:.2f}\")"
        )
    new_header = header.copy()
    new_header.update(w_new.to_header(relax=True))
    new_header["ASTR_SRC"] = ("PS1-ASTROPY", "Astrometric solution source")
    new_header["ASTR_NM"] = (int(np.count_nonzero(good)), "Number of matched sources")
    new_header["ASTR_DRA"] = (float(dra * 3600.0), "Initial median dRA in arcsec")
    new_header["ASTR_DDE"] = (float(ddec * 3600.0), "Initial median dDec in arcsec")
    new_header["ASTR_RMS"] = (fit_rms, "RMS sky residual of final fit (arcsec)")

    if output is None:
        output = os.path.join(os.path.dirname(filename), "a" + os.path.basename(filename))
    fits.writeto(output, data, header=new_header, overwrite=True)
    hdul.close()

    return {
        "file": filename,
        "output": output,
        "n_detected": int(len(x)),
        "n_catalog": int(len(cat_world)),
        "n_coarse": int(np.count_nonzero(coarse_keep)),
        "n_matched": int(np.count_nonzero(good)),
        "n_inlier": int(len(xg)),
        "fit_rms_arcsec": float(fit_rms),
        "used_fine_arcsec": float(used_fine_arcsec),
        "coarse_dra_arcsec": float(dra * 3600.0),
        "coarse_ddec_arcsec": float(ddec * 3600.0),
    }


def parse_args():
    p = argparse.ArgumentParser(description="Solve astrometry with astropy + PS1 catalog.")
    p.add_argument("files", nargs="+", help="Input FITS files.")
    p.add_argument("-o", "--output", default="", help="Output FITS path (single input only).")
    p.add_argument("--band", default="r", help="PS1 band for magnitude filtering (g/r/i/z/y).")
    p.add_argument("--mag-limit", type=float, default=21.5, help="PS1 faint magnitude limit.")
    p.add_argument("--threshold-sigma", type=float, default=5.0, help="DAOStarFinder threshold in sigma.")
    p.add_argument("--coarse-arcsec", type=float, default=120.0, help="Coarse matching radius in arcsec.")
    p.add_argument("--fine-arcsec", type=float, default=2.2, help="Fine matching radius in arcsec.")
    p.add_argument("--min-matches", type=int, default=12, help="Minimum matches required for WCS fit.")
    p.add_argument("--saturation", type=float, default=None, help="Saturation cut for source peaks.")
    return p.parse_args()


def main():
    args = parse_args()
    if args.output and len(args.files) != 1:
        print("Error: --output can only be used with a single input file.")
        return 1

    failed = []
    for f in args.files:
        out = args.output if args.output else None
        try:
            info = solve_one(
                f,
                output=out,
                band=args.band,
                mag_limit=args.mag_limit,
                threshold_sigma=args.threshold_sigma,
                coarse_arcsec=args.coarse_arcsec,
                fine_arcsec=args.fine_arcsec,
                min_matches=args.min_matches,
                saturation=args.saturation,
            )
            print(
                f"Solved {info['file']} -> {info['output']} | "
                f"det={info['n_detected']} cat={info['n_catalog']} "
                f"coarse={info['n_coarse']} matched={info['n_matched']} "
                f"dra={info['coarse_dra_arcsec']:.1f}\" ddec={info['coarse_ddec_arcsec']:.1f}\""
            )
        except Exception as e:
            print(f"Failed {f}: {e}")
            failed.append(f)

    if failed:
        print("\nFailed files:")
        for f in failed:
            print(f"  {f}")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
