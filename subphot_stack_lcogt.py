#!/usr/bin/env python3
"""
subphot_stack_lcogt.py
----------------------
Stack LCOGT (BANZAI e91) exposures into one coadd per date-folder per filter
using SWarp, ready for subphot_subtract.py.

Expected layout (as delivered by the LCOGT archive):

    LCOGT/
      03042026/           <- date folder (any name)
        r/                <- filter sub-folder: g, r, i, u, zs, ...
          tfn1m001-fa20-20260402-0416-e91.fits.fz
          ...
      14032026/
        g/ ...

For every <date>/<filter>/ the script:
  1. flattens the fpack-compressed multi-extension .fz files,
  2. optionally rejects frames with poor seeing (L1FWHM),
  3. coadds them with SWarp centred on the target (CAT-RA/CAT-DEC),
  4. writes <date>/<OBJECT>_<filter>_<YYYYMMDD>_stack.fits with corrected
     EXPTIME/MJD-OBS/GAIN/NCOMBINE keywords.

Usage
-----
    python subphot_stack_lcogt.py -d data/ZTF26aakjzdt_FINAL/LCOGT
    python subphot_stack_lcogt.py -d data/ZTF26aakjzdt_FINAL/LCOGT/03042026 -fb r
    python subphot_stack_lcogt.py -d .../LCOGT --combine MEDIAN --seeing-max 3.0
    python subphot_stack_lcogt.py -d .../LCOGT --dry-run
"""
import argparse, glob, os, re, shutil, sys
import numpy as np
from astropy.io import fits
from astropy.time import Time

from subphot_credentials import path, swarp_path

info_g, warn_y, warn_r = '[STACK]   ::', '[WARNING] ::', '[ERROR]   ::'

FILTER_DIRS = ['u', 'g', 'r', 'i', 'z', 'zs', 'up', 'gp', 'rp', 'ip']
COPY_KEYWORDS = ('OBJECT,TELESCOP,INSTRUME,FILTER,ORIGIN,SITEID,TELID,PROPID,'
                 'CAT-RA,CAT-DEC,RA,DEC,AIRMASS,UTSTART,PIXSCALE,GAIN,RDNOISE,SATURATE')


def flatten_fz(path_in, out_dir):
    """Write a single-HDU copy of a compressed/multi-extension FITS; return its path."""
    hdul = fits.open(path_in)
    try:
        if hdul[0].data is not None:
            return path_in
        for hdu in hdul[1:]:
            if getattr(hdu, 'is_image', False) and hdu.data is not None:
                base = os.path.basename(path_in)
                for suf in ('.fits.fz', '.fits.gz', '.fz', '.fits'):
                    if base.endswith(suf):
                        base = base[:-len(suf)]
                        break
                out = os.path.join(out_dir, base + '_flat.fits')
                fits.PrimaryHDU(data=hdu.data, header=hdu.header.copy()).writeto(out, overwrite=True)
                return out
        raise RuntimeError('no image HDU found')
    finally:
        hdul.close()


def stack_filter_dir(filt_dir, combine, seeing_max, dry_run, keep_tmp):
    date_dir = os.path.dirname(os.path.abspath(filt_dir))
    filt = os.path.basename(os.path.normpath(filt_dir))

    files = sorted(f for f in glob.glob(os.path.join(filt_dir, '*.fits*'))
                   if '_stack' not in f and '_flat' not in f and not f.endswith('.weight.fits'))
    if len(files) == 0:
        print(warn_y + f' {filt_dir}: no FITS files, skipping')
        return None
    print(info_g + f' {filt_dir}: {len(files)} frames')

    tmp_dir = os.path.join(filt_dir, '_stacktmp')
    os.makedirs(tmp_dir, exist_ok=True)

    flats, mjds, exptimes, gains, fwhms = [], [], [], [], []
    hdr0 = None
    try:
        for f in files:
            flat = flatten_fz(f, tmp_dir)
            h = fits.open(flat)[0].header
            fwhm = h.get('L1FWHM', None)
            if seeing_max is not None and fwhm is not None and float(fwhm) > seeing_max:
                print(warn_y + f'   rejecting {os.path.basename(f)}: L1FWHM={float(fwhm):.2f}" > {seeing_max}"')
                continue
            flats.append(flat)
            if hdr0 is None:
                hdr0 = h
            mjds.append(float(h.get('MJD-OBS', h.get('MJD', np.nan))))
            exptimes.append(float(h.get('EXPTIME', 0)))
            gains.append(float(h.get('GAIN', 1.0)))
            if fwhm is not None:
                fwhms.append(float(fwhm))

        if len(flats) == 0:
            print(warn_r + f' {filt_dir}: all frames rejected, skipping')
            return None

        obj = re.sub(r'\s+', '', str(hdr0.get('OBJECT', 'unknown')))
        utdate = re.sub('-', '', str(hdr0.get('DATE-OBS', 'T')).split('T')[0])
        out_name = os.path.join(date_dir, f'{obj}_{filt}_{utdate}_stack.fits')
        weight_name = out_name.replace('.fits', '.weight.fits')

        # centre on the target so all dithers contribute symmetrically
        ra = str(hdr0.get('CAT-RA', hdr0.get('RA', ''))).strip()
        dec = str(hdr0.get('CAT-DEC', hdr0.get('DEC', ''))).strip()
        ps = float(hdr0.get('PIXSCALE', 0.389))
        size = max(int(hdr0.get('NAXIS1', 4096)), int(hdr0.get('NAXIS2', 4096)))

        cmd = (f"{swarp_path} {' '.join(flats)} -c {path}config_files/config_comb.swarp"
               f" -CENTER '{ra} {dec}' -PIXEL_SCALE {ps} -IMAGE_SIZE '{size},{size}'"
               f" -COMBINE_TYPE {combine} -SUBTRACT_BACK Y"
               f" -COPY_KEYWORDS '{COPY_KEYWORDS}'"
               f" -IMAGEOUT_NAME {out_name} -WEIGHTOUT_NAME {weight_name}"
               f" -RESAMPLE_DIR {tmp_dir} -VERBOSE_TYPE QUIET")
        if dry_run:
            print(info_g + f'   would stack {len(flats)} frames -> {out_name}')
            return None

        status = os.system(cmd)
        if status != 0 or not os.path.exists(out_name):
            print(warn_r + f' SWarp failed (status {status}) for {filt_dir}')
            return None

        # fix up the coadd header for downstream photometry
        with fits.open(out_name, mode='update') as hdul:
            h = hdul[0].header
            n = len(flats)
            h['NCOMBINE'] = (n, 'number of frames combined')
            h['MJD-OBS'] = (float(np.nanmean(mjds)), 'mean MJD of combined frames')
            h['DATE-OBS'] = (Time(float(np.nanmean(mjds)), format='mjd').isot,
                             'UT date at mean MJD')
            # AVERAGE coadd keeps single-frame count levels: EXPTIME stays the
            # mean single exposure; the summed time is recorded separately and
            # the effective gain scales with the number of averaged frames
            h['EXPTIME'] = (float(np.nanmean(exptimes)), 'mean single-frame exposure [s]')
            h['TOTEXPT'] = (float(np.nansum(exptimes)), 'summed exposure time [s]')
            # belt-and-braces: header_kw needs these; fall back to the first input if
            # SWarp's COPY_KEYWORDS ever misses one
            for kw in ('UTSTART', 'CAT-RA', 'CAT-DEC', 'AIRMASS', 'PIXSCALE', 'FILTER'):
                if kw not in h and kw in hdr0:
                    h[kw] = hdr0[kw]
            if combine.upper() == 'AVERAGE':
                h['GAIN'] = (float(np.nanmean(gains)) * n, f'effective gain (AVERAGE of {n})')
            if len(fwhms) > 0:
                h['L1FWHM'] = (float(np.nanmedian(fwhms)), 'median L1FWHM of inputs [arcsec]')
        print(info_g + f'   stacked {len(flats)} frames -> {out_name}'
              + (f'  (median L1FWHM {np.nanmedian(fwhms):.2f}")' if fwhms else ''))
        return out_name
    finally:
        if not keep_tmp:
            shutil.rmtree(tmp_dir, ignore_errors=True)


def main():
    p = argparse.ArgumentParser(description='Stack LCOGT exposures per date folder per filter with SWarp')
    p.add_argument('-d', '--dirs', nargs='+', required=True,
                   help='LCOGT top-level directory (containing date folders) or individual date folder(s)')
    p.add_argument('-fb', '--filters', nargs='+', default=None,
                   help='only stack these filter sub-folders (default: all found)')
    p.add_argument('--combine', default='AVERAGE', choices=['AVERAGE', 'MEDIAN', 'CLIPPED', 'WEIGHTED'],
                   help='SWarp COMBINE_TYPE (default AVERAGE)')
    p.add_argument('--seeing-max', type=float, default=None,
                   help='reject frames with L1FWHM above this value in arcsec (default: keep all)')
    p.add_argument('--dry-run', action='store_true', help='list what would be stacked without running SWarp')
    p.add_argument('--keep-tmp', action='store_true', help='keep the flattened temporary files')
    args = p.parse_args()

    # resolve to a list of filter sub-folders
    filter_dirs = []
    for d in args.dirs:
        if not os.path.isdir(d):
            print(warn_r + f' {d} is not a directory'); continue
        subs = sorted(os.listdir(d))
        if any(os.path.basename(s) in FILTER_DIRS for s in subs):        # d is a date folder
            date_dirs = [d]
        else:                                                            # d is the LCOGT top level
            date_dirs = [os.path.join(d, s) for s in subs if os.path.isdir(os.path.join(d, s))]
        for dd in date_dirs:
            for s in sorted(os.listdir(dd)):
                fd = os.path.join(dd, s)
                if os.path.isdir(fd) and s in FILTER_DIRS:
                    if args.filters is None or s in args.filters or s.rstrip('ps') in (args.filters or []):
                        filter_dirs.append(fd)

    if len(filter_dirs) == 0:
        print(warn_r + ' no filter sub-folders found'); sys.exit(1)

    print(info_g + f' {len(filter_dirs)} filter folder(s) to stack')
    stacks = []
    for fd in filter_dirs:
        out = stack_filter_dir(fd, args.combine, args.seeing_max, args.dry_run, args.keep_tmp)
        if out: stacks.append(out)
    print(info_g + f' done: {len(stacks)} stack(s) written')
    for s in stacks: print(info_g + f'   {s}')


if __name__ == '__main__':
    main()
