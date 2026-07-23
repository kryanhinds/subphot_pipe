#!/usr/bin/env python3
"""
stack_by_night.py
-----------------
Group FITS images in a directory by night and filter, then co-add each group
using astroalign (star-pattern matching, no WCS required) + numpy mean stacking.

Usage
-----
    # Stack everything in a folder
    python stack_by_night.py -f data/ZTF26aakjzdt/

    # Stack only images from nights after the 24th of any month
    python stack_by_night.py -f data/ZTF26aakjzdt/ --after-day 24

    # Dry run — show groups without stacking
    python stack_by_night.py -f data/ZTF26aakjzdt/ --dry-run

    # Use median instead of mean combine
    python stack_by_night.py -f data/ZTF26aakjzdt/ --combine median

Notes
-----
* astroalign aligns by matching triangle patterns of bright stars — the WCS
  does not need to be accurate or even present.
* Groups on (astro-night, OBJECT, FILTER).  Astro-night = UT date − 12 h so
  images taken before midnight UT are assigned to the previous calendar night.
* Requires ≥ 2 images per group (single images are skipped).
* Output filename:  stacked/<OBJECT>_<NIGHT>_<FILTER>_stacked.fits
* The output header is copied from the first (reference) image with EXPTIME
  set to the total summed exposure time and NCOMBINE added.
"""

import argparse
import os
import re
import shutil
import tempfile
from collections import defaultdict
from datetime import datetime, timedelta

import numpy as np
from astropy.io import fits

# ── optional credentials ──────────────────────────────────────────────────────
try:
    from subphot_credentials import path as _cred_path
except ImportError:
    _cred_path = os.path.dirname(os.path.abspath(__file__)) + '/'


# ─────────────────────────────────────────────────────────────────────────────
def parse_date_obs(s):
    s = str(s).strip()
    for fmt in ('%Y-%m-%dT%H:%M:%S.%f', '%Y-%m-%dT%H:%M:%S',
                '%Y-%m-%d %H:%M:%S', '%Y-%m-%d'):
        try:
            return datetime.strptime(s[:len(fmt) + 4], fmt)
        except ValueError:
            continue
    raise ValueError(f'Cannot parse date string: {s!r}')


def astro_night(dt):
    """UT date at noon: observations in the same night share the same label."""
    return (dt - timedelta(hours=12)).strftime('%Y-%m-%d')


def read_kw(hdr, *keys, default=None):
    for k in keys:
        if k in hdr:
            return str(hdr[k]).strip()
    return default


# ─────────────────────────────────────────────────────────────────────────────
def collect_images(folder, after_day=None):
    """Scan *folder* for FITS files; return dict keyed by (night, obj, filt)."""
    groups = defaultdict(list)
    fnames = sorted(
        f for f in os.listdir(folder)
        if f.lower().endswith(('.fits', '.fit'))
        and not f.startswith('.')
        and not f.endswith('tpv.fits')
    )
    if not fnames:
        print(f'[stack] No FITS files found in {folder}')
        return groups

    for fname in fnames:
        fpath = os.path.join(folder, fname)
        try:
            with fits.open(fpath, memmap=True) as hdul:
                hdr = hdul[0].header
                if hdr.get('NAXIS', 0) < 2:
                    continue

                date_str = read_kw(hdr, 'DATE-OBS', 'DATE', 'UTC', 'MJD_OBS')
                if date_str is None:
                    print(f'  [skip] {fname}: no date keyword')
                    continue

                # MJD as float
                if re.match(r'^\d{5}', date_str):
                    from astropy.time import Time as ATime
                    dt = ATime(float(date_str), format='mjd').to_datetime()
                else:
                    dt = parse_date_obs(date_str)

                night = astro_night(dt)

                if after_day is not None:
                    if int(night.split('-')[2]) <= after_day:
                        continue

                filt = read_kw(hdr, 'FILTER1', 'FILTER', 'NCFLTNM2',
                               'SEQID', default='unknown')
                filt = re.sub(r'[-_\s].*', '', filt).strip()

                obj = read_kw(hdr, 'OBJECT', 'TARGET', 'TCSTGT',
                              default='unknown')
                obj = re.sub(r'[\s/]', '_', obj).strip()

                groups[(night, obj, filt)].append(fpath)

        except Exception as exc:
            print(f'  [skip] {fname}: {exc}')

    return groups


# ─────────────────────────────────────────────────────────────────────────────
def _align_astroalign(src, ref):
    """Try astroalign with progressively relaxed detection thresholds.

    Tries detection_sigma = 3, 2, 1 in turn, increasing max_control_points
    each time to maximise the chance of finding enough star triangles on
    sparse fields (SEDM-P60 images often have only 3–5 bright sources).
    Returns the registered array, or raises if all attempts fail.
    """
    import astroalign as aa

    attempts = [
        dict(detection_sigma=3, max_control_points=50,  min_area=5),
        dict(detection_sigma=2, max_control_points=100, min_area=3),
        dict(detection_sigma=1, max_control_points=200, min_area=2),
    ]
    last_exc = None
    for kwargs in attempts:
        try:
            registered, _ = aa.register(src, ref, **kwargs)
            return registered, f'astroalign (sigma={kwargs["detection_sigma"]})'
        except Exception as exc:
            last_exc = exc
    raise last_exc


def _align_crosscorr(src, ref):
    """Shift-only alignment via phase cross-correlation.

    Used as a fallback when astroalign cannot find enough star triangles.
    Handles only integer-pixel shifts — sufficient for same-telescope images
    taken on different nights where pointing differences are small.
    """
    from scipy.ndimage import shift as nd_shift
    try:
        from skimage.registration import phase_cross_correlation
        shifts, _, _ = phase_cross_correlation(ref, src, upsample_factor=10)
    except ImportError:
        # older scikit-image API
        from skimage.feature import register_translation
        shifts, _, _ = register_translation(ref, src, upsample_factor=10)

    registered = nd_shift(src, shifts)
    return registered, f'phase_cross_correlation (shift={shifts[0]:.1f},{shifts[1]:.1f}px)'


def align_and_stack(image_paths, combine='mean'):
    """Align all images to the first, then combine.

    Alignment strategy (per image):
      1. astroalign with detection_sigma 3→2→1 (triangle pattern matching,
         no WCS required, handles rotation/scale differences)
      2. Fallback: phase cross-correlation shift (for very sparse fields where
         astroalign cannot find enough star triangles)

    Parameters
    ----------
    image_paths : list of str
    combine     : 'mean' or 'median'

    Returns
    -------
    stacked_data  : 2-D numpy array
    ref_header    : FITS header from the reference (first) image
    total_exptime : float
    n_combined    : int
    """
    ref_hdul   = fits.open(image_paths[0])
    ref_data   = ref_hdul[0].data.astype(np.float64)
    ref_header = ref_hdul[0].header.copy()
    ref_hdul.close()

    aligned_frames = [ref_data]
    total_exptime  = float(ref_header.get('EXPTIME', 0))
    failed         = []

    for fpath in image_paths[1:]:
        fname = os.path.basename(fpath)
        try:
            with fits.open(fpath) as hdul:
                src_data = hdul[0].data.astype(np.float64)
                exptime  = float(hdul[0].header.get('EXPTIME', 0))

            # ── attempt 1: astroalign ──────────────────────────────────────
            try:
                registered, method = _align_astroalign(src_data, ref_data)
            except Exception as aa_exc:
                # ── attempt 2: phase cross-correlation ────────────────────
                try:
                    registered, method = _align_crosscorr(src_data, ref_data)
                    print(f'    [align warn] {fname}: astroalign failed '
                          f'({aa_exc}), used fallback {method}')
                except Exception as cc_exc:
                    raise RuntimeError(
                        f'astroalign: {aa_exc}  |  crosscorr: {cc_exc}'
                    ) from None

            else:
                print(f'    [align ok  ] {fname}  [{method}]')

            aligned_frames.append(registered)
            total_exptime += exptime

        except Exception as exc:
            print(f'    [align FAIL] {fname}: {exc}')
            failed.append(fname)

    if len(aligned_frames) < 2:
        raise RuntimeError(
            f'Only {len(aligned_frames)} frame(s) aligned successfully '
            f'(need ≥ 2).  Failed: {failed}'
        )

    stack = np.array(aligned_frames)
    stacked_data = np.median(stack, axis=0) if combine == 'median' \
                   else np.mean(stack, axis=0)

    return stacked_data, ref_header, total_exptime, len(aligned_frames)


# ─────────────────────────────────────────────────────────────────────────────
def write_stacked(out_path, data, ref_header, total_exptime, n_combined):
    """Write stacked array to *out_path* with a cleaned-up header."""
    hdr = ref_header.copy()

    # Remove SIP distortion (no longer valid after pixel-space alignment)
    for key in list(hdr.keys()):
        if re.match(r'^[AB]P?_\d', key) or key in ('A_ORDER', 'B_ORDER',
                                                     'AP_ORDER', 'BP_ORDER'):
            del hdr[key]
    for ct in ('CTYPE1', 'CTYPE2'):
        if ct in hdr:
            hdr[ct] = re.sub(r'-SIP$', '', hdr[ct])

    hdr['EXPTIME']  = (total_exptime,  'Total stacked exposure time [s]')
    hdr['NCOMBINE'] = (n_combined,     'Number of frames co-added')
    hdr['STCKTYPE'] = ('astroalign',   'Stacking method')

    out_dir = os.path.dirname(out_path)
    tmp_fd, tmp_path = tempfile.mkstemp(dir=out_dir, suffix='.fits')
    os.close(tmp_fd)
    try:
        fits.writeto(tmp_path, data.astype(np.float32), hdr, overwrite=True)
        shutil.move(tmp_path, out_path)
    except Exception:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise


# ─────────────────────────────────────────────────────────────────────────────
def stack_folder(folder, after_day=None, dry_run=False,
                 min_images=2, combine='mean'):
    folder      = os.path.abspath(folder)
    stacked_dir = os.path.join(folder, 'stacked')

    if not dry_run:
        os.makedirs(stacked_dir, exist_ok=True)

    groups = collect_images(folder, after_day=after_day)

    if not groups:
        print('[stack] No groups to stack (after filtering).')
        return

    print(f'\n[stack] Found {len(groups)} group(s) in {folder}')
    print(f'[stack] Output directory : {stacked_dir}')
    print(f'[stack] Combine method   : {combine}\n')

    n_ok = n_skip = 0

    for (night, obj, filt), images in sorted(groups.items()):
        if len(images) < min_images:
            print(f'  [skip] {night} / {obj} / {filt}: '
                  f'{len(images)} image (need ≥ {min_images})')
            n_skip += 1
            continue

        out_name = f'{obj}_{night}_{filt}_stacked.fits'
        out_path = os.path.join(stacked_dir, out_name)

        print(f'  [{night}]  {obj}  filter={filt}  ({len(images)} images)')
        for im in images:
            print(f'    + {os.path.basename(im)}')
        print(f'    → {out_name}')

        if dry_run:
            continue

        try:
            data, hdr, exptime, n = align_and_stack(images, combine=combine)
            write_stacked(out_path, data, hdr, exptime, n)
            print(f'  [ok]  Wrote {out_path}  '
                  f'(EXPTIME={exptime:.0f}s, N={n})\n')
            n_ok += 1
        except Exception as exc:
            print(f'  [FAIL] {out_name}: {exc}\n')
            n_skip += 1

    print(f'[stack] Done — {n_ok} stack(s) written, {n_skip} skipped/failed.')


# ─────────────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    ap = argparse.ArgumentParser(
        description='Align with astroalign and co-add FITS images by night.')
    ap.add_argument('-f', '--folder', required=True,
                    help='Directory containing science FITS images')
    ap.add_argument('--after-day', type=int, default=None, metavar='DAY',
                    help='Only stack images from nights after this day-of-month '
                         '(e.g. 24 → keep nights 25th onwards)')
    ap.add_argument('--dry-run', action='store_true',
                    help='Show groups without stacking')
    ap.add_argument('--min-images', type=int, default=2,
                    help='Minimum images per group (default 2)')
    ap.add_argument('--combine', choices=['mean', 'median'], default='mean',
                    help='Combine method (default: mean)')
    args = ap.parse_args()

    stack_folder(
        folder     = args.folder,
        after_day  = args.after_day,
        dry_run    = args.dry_run,
        min_images = args.min_images,
        combine    = args.combine,
    )
