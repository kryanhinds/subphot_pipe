#!/usr/bin/env python3
"""
Create color composite images from multi-filter FITS files.

Usage:
    python make_color_image.py file1.fits file2.fits file3.fits \
        -wl g:473 r:622 i:763 z:905 -o output.png

The script automatically maps filters to RGB based on wavelengths (red=longest, green=middle, blue=shortest).
Supports various image scaling methods for optimal contrast.
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, Tuple, List, Optional

import numpy as np
from astropy.io import fits
from astropy.visualization import ImageNormalize, AsinhStretch, PercentileInterval, LinearStretch
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from PIL import Image


def read_fits_image(filepath: str) -> Tuple[np.ndarray, dict]:
    """
    Read FITS file and extract image data and header.

    Parameters
    ----------
    filepath : str
        Path to FITS file

    Returns
    -------
    image : ndarray
        2D image data
    header : dict
        FITS header as dictionary
    """
    with fits.open(filepath) as hdul:
        image = hdul[0].data
        header = dict(hdul[0].header)

    return image, header


def extract_filter(header: dict, filter_keyword: str = 'FILTER') -> str:
    """
    Extract filter name from FITS header.

    Parameters
    ----------
    header : dict
        FITS header dictionary
    filter_keyword : str
        FITS keyword to search for filter (default: 'FILTER')

    Returns
    -------
    filter_name : str
        Single character filter name (g, r, i, z, etc.)
    """
    # Try common filter keywords
    keywords = [filter_keyword, 'FILTER', 'FILT', 'FILTID', 'FILTERID', 'INSFLTR']

    for kw in keywords:
        if kw in header:
            filt = str(header[kw]).strip().lower()
            # Extract single character if it's embedded in a longer string
            if len(filt) > 1:
                for char in filt:
                    if char in 'grizuy':
                        return char
            elif filt in 'grizuy':
                return filt

    raise ValueError(f"Could not identify filter in header. Tried keywords: {keywords}")


def scale_image(image: np.ndarray, method: str = 'asinh') -> np.ndarray:
    """
    Scale image data to 0-255 range for visualization.

    Parameters
    ----------
    image : ndarray
        Raw image data
    method : str
        Scaling method: 'asinh', 'percentile', 'linear', 'sqrt'

    Returns
    -------
    scaled : ndarray
        Scaled image (0-255)
    """
    # Remove NaNs and Infs for scaling
    valid_mask = np.isfinite(image)
    if not np.any(valid_mask):
        raise ValueError("Image contains no valid (finite) data")

    vmin, vmax = np.nanpercentile(image[valid_mask], [1, 99])

    if method == 'asinh':
        # Asinh stretch - good for showing both faint and bright details
        norm = ImageNormalize(image, interval=PercentileInterval(99.5),
                            stretch=AsinhStretch(a=0.02))
    elif method == 'percentile':
        # Simple percentile stretch
        norm = Normalize(vmin=vmin, vmax=vmax, clip=True)
    elif method == 'sqrt':
        # Square root stretch
        image_sqrt = np.sqrt(np.clip(image, vmin, vmax) - vmin)
        return (image_sqrt / np.nanmax(image_sqrt) * 255).astype(np.uint8)
    elif method == 'linear':
        # Linear stretch
        norm = Normalize(vmin=vmin, vmax=vmax, clip=True)
    else:
        raise ValueError(f"Unknown scaling method: {method}")

    scaled = norm(image)
    scaled = (np.clip(scaled, 0, 1) * 255).astype(np.uint8)

    return scaled


def align_images(images: Dict[str, np.ndarray],
                 reference_filter: Optional[str] = None) -> Dict[str, np.ndarray]:
    """
    Align images to a common grid (assumes WCS headers are present in FITS).

    For now, assumes images are already aligned or will check if they're the same shape.
    For full alignment, would need reproject_interp from reproject module.

    Parameters
    ----------
    images : dict
        Filter -> image array mapping
    reference_filter : str, optional
        Which filter to use as reference. If None, uses first one.

    Returns
    -------
    aligned : dict
        Aligned images
    """
    shapes = [img.shape for img in images.values()]

    if len(set(shapes)) > 1:
        print("Warning: Images have different shapes:")
        for filt, shape in zip(images.keys(), shapes):
            print(f"  {filt}: {shape}")
        print("\nAssuming images are pre-aligned with WCS.")
        print("For automatic alignment, ensure FITS headers contain valid WCS.")

    return images


def map_filters_to_rgb(filters: List[str],
                       wavelengths: Dict[str, float]) -> Dict[str, str]:
    """
    Intelligently map filters to RGB channels based on wavelengths.

    Parameters
    ----------
    filters : list
        List of available filter names (e.g., ['g', 'r', 'i'])
    wavelengths : dict
        Filter -> wavelength (nm) mapping

    Returns
    -------
    mapping : dict
        Filter -> 'R'/'G'/'B' channel mapping
    """
    if len(filters) < 3:
        raise ValueError(f"Need at least 3 filters for RGB, got {len(filters)}: {filters}")

    # Sort filters by wavelength
    sorted_filters = sorted(filters, key=lambda f: wavelengths[f])

    mapping = {
        sorted_filters[-1]: 'R',  # Longest wavelength -> Red
        sorted_filters[1]: 'G',   # Middle wavelength -> Green
        sorted_filters[0]: 'B',   # Shortest wavelength -> Blue
    }

    print(f"Filter -> RGB mapping (by wavelength):")
    for filt, channel in sorted(mapping.items(), key=lambda x: wavelengths[x[0]]):
        print(f"  {filt} ({wavelengths[filt]} nm) -> {channel}")

    return mapping


def create_color_image(files: List[str],
                      wavelengths: Dict[str, float],
                      output: str = 'color_image.png',
                      scaling_method: str = 'asinh',
                      filter_keyword: str = 'FILTER') -> None:
    """
    Create a color composite from multi-filter FITS files.

    Parameters
    ----------
    files : list
        List of FITS file paths
    wavelengths : dict
        Filter -> wavelength (nm) mapping for RGB assignment
    output : str
        Output filename
    scaling_method : str
        Image scaling method ('asinh', 'percentile', 'linear', 'sqrt')
    filter_keyword : str
        FITS header keyword for filter name
    """

    print(f"Processing {len(files)} files...")

    # Read and catalog images
    images_by_filter = {}
    headers_by_filter = {}

    for filepath in files:
        print(f"  Reading {Path(filepath).name}...", end=' ')
        image, header = read_fits_image(filepath)
        filt = extract_filter(header, filter_keyword)

        if filt in images_by_filter:
            print(f"WARNING: {filt} already loaded, skipping duplicate")
            continue

        images_by_filter[filt] = image
        headers_by_filter[filt] = header
        print(f"[{filt}] {image.shape}")

    # Verify we have the filters in the wavelength dict
    available_filters = list(images_by_filter.keys())
    for filt in available_filters:
        if filt not in wavelengths:
            raise ValueError(f"Filter '{filt}' found in images but not in wavelengths dict")

    print(f"\nAvailable filters: {available_filters}")

    # Align images (currently just checks shapes)
    images_by_filter = align_images(images_by_filter)

    # Map filters to RGB
    rgb_mapping = map_filters_to_rgb(available_filters, wavelengths)

    # Create RGB array
    ref_shape = list(images_by_filter.values())[0].shape
    rgb_image = np.zeros((*ref_shape, 3), dtype=np.uint8)

    channel_map = {'R': 0, 'G': 1, 'B': 2}

    print(f"\nScaling images using '{scaling_method}' method...")
    for filt, image in images_by_filter.items():
        channel = rgb_mapping[filt]
        channel_idx = channel_map[channel]

        scaled = scale_image(image, method=scaling_method)
        rgb_image[:, :, channel_idx] = scaled
        print(f"  {filt} -> {channel}: range [{scaled.min()}, {scaled.max()}]")

    # Save output
    print(f"\nSaving to {output}...")
    img = Image.fromarray(rgb_image, mode='RGB')
    img.save(output)
    print(f"✓ Color image saved: {output}")
    print(f"  Size: {img.size}")

    return rgb_image


def main():
    parser = argparse.ArgumentParser(
        description='Create color composite images from multi-filter FITS files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Combine g, r, i images with default wavelengths
  python make_color_image.py image_g.fits image_r.fits image_i.fits -o color.png

  # Specify custom wavelengths (in nm)
  python make_color_image.py g.fits r.fits i.fits \\
    -wl g:473 r:622 i:763 -o custom.png

  # Include z-band (will use g, r, i as RGB)
  python make_color_image.py g.fits r.fits i.fits z.fits \\
    -wl g:473 r:622 i:763 z:905 -o rgbi.png

  # Use percentile scaling instead of asinh
  python make_color_image.py *.fits -wl g:473 r:622 i:763 -s percentile
        """)

    parser.add_argument('files', nargs='+', help='FITS files to combine')
    parser.add_argument('-o', '--output', default='color_image.png',
                       help='Output filename (default: color_image.png)')
    parser.add_argument('-wl', '--wavelengths', nargs='+',
                       help='Filter wavelengths as "filter:wavelength_nm" '
                            '(default: g:473 r:622 i:763 z:905)')
    parser.add_argument('-s', '--scaling', default='asinh',
                       choices=['asinh', 'percentile', 'linear', 'sqrt'],
                       help='Image scaling method (default: asinh)')
    parser.add_argument('-fk', '--filter-keyword', default='FILTER',
                       help='FITS header keyword for filter (default: FILTER)')

    args = parser.parse_args()

    # Default wavelengths (SDSS/ZTF bandpasses, central wavelengths in nm)
    default_wavelengths = {
        'g': 473,
        'r': 622,
        'i': 763,
        'z': 905,
        'y': 960,
        'u': 355,
    }

    # Parse user wavelengths if provided
    wavelengths = default_wavelengths.copy()
    if args.wavelengths:
        for wl_spec in args.wavelengths:
            try:
                filt, wl = wl_spec.split(':')
                wavelengths[filt.lower()] = float(wl)
            except ValueError:
                print(f"Error parsing wavelength spec: {wl_spec}")
                print("Expected format: 'filter:wavelength_nm' (e.g., 'g:473')")
                sys.exit(1)

    print("=" * 60)
    print("Color Image Creator for Astronomical FITS Files")
    print("=" * 60)
    print(f"Wavelengths (nm): {wavelengths}\n")

    try:
        create_color_image(
            files=args.files,
            wavelengths=wavelengths,
            output=args.output,
            scaling_method=args.scaling,
            filter_keyword=args.filter_keyword
        )
    except Exception as e:
        print(f"\nError: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
