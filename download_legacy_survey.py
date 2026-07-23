#!/usr/bin/env python3
"""
Download Legacy Survey reference images for subphot_pipe

Usage:
    python download_legacy_survey.py 219.317292 71.841750 -b g r z -o ref_imgs/
"""

import requests
import gzip
import shutil
import argparse
from pathlib import Path
from astropy.io import fits

def download_legacy_survey(ra, dec, bands='g', output_dir='ref_imgs',
                          size=800, pixscale=0.262, layer='ls-dr9'):
    """
    Download Legacy Survey cutout image

    Parameters:
    -----------
    ra : float
        Right ascension in degrees
    dec : float
        Declination in degrees
    bands : str
        Filter bands: 'g', 'r', 'i', 'z' (can combine: 'gr', 'gri', etc.)
    output_dir : str
        Directory to save FITS file
    size : int
        Image size in pixels (default 800 = 400x400 arcmin at 0.262 arcsec/pix)
    pixscale : float
        Pixel scale in arcsec/pixel (default 0.262)
    layer : str
        Legacy Survey data release (default 'ls-dr9')

    Returns:
    --------
    str : Path to FITS file, or None if failed
    """

    # Ensure output directory exists
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # URL for Legacy Survey FITS cutout using viewer API
    # This is the correct endpoint from the Legacy Survey documentation
    url = (f"https://www.legacysurvey.org/viewer/cutout.fits?"
           f"ra={ra:.6f}&dec={dec:.6f}&"
           f"layer={layer}&"
           f"bands={bands}&"
           f"size={size}&"
           f"pixscale={pixscale}")

    output_name = f"legacysurvey_{''.join(bands)}.fits"
    output_path = Path(output_dir) / output_name

    print(f"\n{'='*70}")
    print(f"Downloading {','.join(list(bands))}-band Legacy Survey ({size}x{size} px)")
    print(f"{'='*70}")
    print(f"  RA={ra:.6f}, Dec={dec:.6f}")
    print(f"  Pixel scale: {pixscale} arcsec/px")
    print(f"  Output: {output_path}")

    try:
        # Download
        print("  Downloading...", end='', flush=True)
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        print(" ✓")

        # Check if response is gzipped
        if response.headers.get('content-encoding') == 'gzip' or url.endswith('.gz'):
            print("  Uncompressing...", end='', flush=True)
            import io
            with gzip.open(io.BytesIO(response.content), 'rb') as f:
                with open(output_path, 'wb') as out_f:
                    out_f.write(f.read())
            print(" ✓")
        else:
            # Save directly (already uncompressed)
            with open(output_path, 'wb') as f:
                f.write(response.content)

        # Verify and report
        hdul = fits.open(output_path)
        data = hdul[0].data
        shape = data.shape
        min_val = data.min()
        max_val = data.max()
        mean_val = data.mean()
        print(f"  ✓ Verified: shape {shape}")
        print(f"           range [{min_val:.2f}, {max_val:.2f}]")
        print(f"           mean={mean_val:.2f}")
        hdul.close()

        return str(output_path)

    except requests.exceptions.ConnectionError:
        print(" ✗")
        print(f"  ERROR: Connection failed - check internet or URL")
        return None
    except requests.exceptions.HTTPError as e:
        print(" ✗")
        print(f"  ERROR: {e}")
        print(f"  This image may not exist for these coordinates")
        return None
    except Exception as e:
        print(" ✗")
        print(f"  ERROR: {e}")
        return None

def main():
    parser = argparse.ArgumentParser(
        description='Download Legacy Survey (DR9+) reference images',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download g-band (800x800 pixels)
  python download_legacy_survey.py 219.317292 71.841750 -b g

  # Download g, r bands combined into single FITS
  python download_legacy_survey.py 219.317292 71.841750 -b gr

  # Download larger image (1200x1200 pixels)
  python download_legacy_survey.py 219.317292 71.841750 -b g --size 1200

  # Save to custom directory
  python download_legacy_survey.py 219.317292 71.841750 -b g -o my_refs/
        """
    )

    parser.add_argument('ra', type=float, help='Right ascension (degrees)')
    parser.add_argument('dec', type=float, help='Declination (degrees)')
    parser.add_argument('-b', '--bands', default='g',
                       help='Filter bands: g, r, i, z (combine: gr, gri, etc.) (default: g)')
    parser.add_argument('--size', type=int, default=800,
                       help='Image size in pixels (default: 800 = ~7 arcmin at native scale)')
    parser.add_argument('--pixscale', type=float, default=0.262,
                       help='Pixel scale in arcsec/pixel (default: 0.262)')
    parser.add_argument('-o', '--output', default='ref_imgs',
                       help='Output directory (default: ref_imgs)')
    parser.add_argument('--layer', default='ls-dr9',
                       help='Legacy Survey layer (default: ls-dr9). Options: ls-dr8, ls-dr9, etc.')

    args = parser.parse_args()

    print(f"\n{'='*70}")
    print(f"Legacy Survey Image Downloader")
    print(f"{'='*70}")
    print(f"Target: RA={args.ra:.6f}, Dec={args.dec:.6f}")
    print(f"Bands: {args.bands}")
    print(f"Size: {args.size}x{args.size} px (pixscale={args.pixscale} as/px)")
    print(f"Layer: {args.layer}")
    print(f"Output dir: {args.output}")

    result = download_legacy_survey(args.ra, args.dec, bands=args.bands,
                                    output_dir=args.output, size=args.size,
                                    pixscale=args.pixscale, layer=args.layer)

    print(f"\n{'='*70}")
    if result:
        print(f"✓ Successfully downloaded to: {result}")
        print(f"{'='*70}")
        print(f"\nUse with subphot_subtract.py:")
        print(f"  -refimg {result}")
    else:
        print(f"✗ Download failed")
        print(f"{'='*70}")

if __name__ == '__main__':
    main()
