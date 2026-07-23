# Color Image Creation from Multi-Filter Data

This notebook demonstrates how to create color composite images from a DataFrame of multi-filter FITS files using WCS-based alignment.

## Setup

```python
import pandas as pd
import numpy as np
from color_image_notebook import ColorImageMaker

# Optional: configure matplotlib for notebooks
import matplotlib.pyplot as plt
%matplotlib inline
```

## Example 1: Basic Usage

Suppose you have a DataFrame with your FITS file information:

```python
# Create or load your data
df = pd.DataFrame({
    'filter': ['g', 'r', 'i', 'g', 'r', 'i'],
    'mjd': [59800.123, 59800.125, 59800.127, 59800.128, 59800.130, 59800.132],
    'instrument': ['NOT', 'NOT', 'NOT', 'LT', 'LT', 'LT'],
    'filename': [
        '/path/to/image_NOT_g.fits',
        '/path/to/image_NOT_r.fits',
        '/path/to/image_NOT_i.fits',
        '/path/to/image_LT_g.fits',
        '/path/to/image_LT_r.fits',
        '/path/to/image_LT_i.fits',
    ]
})

# Select best seeing images
df_good = df[df['seeing'] < 1.5]  # or your quality metric

print(df_good)
```

## Example 2: Create Color Image

```python
# Initialize the maker with default wavelengths
maker = ColorImageMaker(
    wavelengths={'g': 473, 'r': 622, 'i': 763},  # Optional custom wavelengths
    scaling_method='asinh',  # or 'percentile', 'linear', 'sqrt'
    verbose=True
)

# Create from DataFrame
rgb_image, info = maker.create_from_dataframe(
    df_good,
    output_path='color_composite.png',
    check_mjd_range=0.01  # Warn if MJD spread > 0.01 days
)

print("Image info:", info)
```

## Example 3: Display in Notebook

```python
# Display the color image
fig = maker.plot(
    figsize=(12, 10),
    title='NOT + LT Color Composite (g, r, i)',
    show_info=True
)
plt.show()
```

## Example 4: Custom Wavelengths

If you have non-standard filters or want to use actual measured wavelengths:

```python
custom_wl = {
    'g': 480,  # Your custom g-band center
    'r': 630,  # Your custom r-band center
    'i': 770,  # Your custom i-band center
}

maker = ColorImageMaker(wavelengths=custom_wl)
rgb_image, info = maker.create_from_dataframe(df_good)
maker.plot()
```

## Example 5: Different Scaling Methods

Try different contrast stretches:

```python
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

for ax, method in zip(axes.flat, ['asinh', 'percentile', 'linear', 'sqrt']):
    maker = ColorImageMaker(scaling_method=method)
    rgb, _ = maker.create_from_dataframe(df_good)
    ax.imshow(rgb)
    ax.set_title(f'{method} scaling')
    ax.axis('off')

plt.tight_layout()
plt.show()
```

## Example 6: Work with Filtered Subsets

Process only certain instruments or filter combinations:

```python
# Only NOT data
df_not = df[df['instrument'] == 'NOT']

# Only certain filters
df_gri = df[df['filter'].isin(['g', 'r', 'i'])]

# Best epochs (within 0.005 days)
mjd_center = df['mjd'].median()
df_epoch = df[(df['mjd'] - mjd_center).abs() < 0.005]

# Combine filters
maker = ColorImageMaker()
rgb, info = maker.create_from_dataframe(df_epoch)
maker.plot()
```

## Example 7: Save Different Formats

```python
# Save as PNG (default, lossless)
maker.save('color_image.png')

# Save as JPEG (lossy, smaller)
from PIL import Image
Image.fromarray(maker.rgb_image).save('color_image.jpg', quality=95)

# Get raw numpy array for further processing
raw_rgb = maker.rgb_image
print(f"RGB shape: {raw_rgb.shape}, dtype: {raw_rgb.dtype}")
```

## Example 8: Batch Process Multiple Objects

```python
objects = ['ZTF25aceeneu', 'ZTF26abc123', 'AT2026xyz']

for obj in objects:
    # Load data for this object
    df_obj = pd.read_csv(f'observations/{obj}_data.csv')
    
    # Keep only good quality images near same epoch
    df_best = df_obj[(df_obj['seeing'] < 1.5) & 
                     ((df_obj['mjd'] - df_obj['mjd'].median()).abs() < 0.01)]
    
    if len(df_best) < 3:
        print(f"Not enough data for {obj}")
        continue
    
    # Create color image
    maker = ColorImageMaker()
    try:
        rgb, info = maker.create_from_dataframe(
            df_best,
            output_path=f'output/{obj}_color.png'
        )
        print(f"✓ Created {obj}")
    except Exception as e:
        print(f"✗ Failed {obj}: {e}")
```

## Notes on WCS Alignment

The script uses **reproject_interp** from the `reproject` package to align all images to a common WCS solution:

1. The **middle filter (by wavelength) is used as reference** for stability
2. All other images are reprojected to match its WCS and shape
3. Missing or invalid data (outside image bounds after reprojection) becomes 0

If alignment fails:
- Check that FITS headers contain valid WCS (`CRPIX`, `CRVAL`, `CD` matrix, etc.)
- Consider using the `sedm_subtract.py` pipeline's astrometry first if WCS is poor
- Or try the simpler command-line version without WCS alignment

## Troubleshooting

**"No valid WCS in header"**
- Your FITS files need astrometric solutions first
- Use `sedm_subtract.py -reastrom` to solve or improve WCS

**"Images have different shapes after alignment"**
- This is usually OK—reproject handles it by padding/cropping
- But very different image sizes may indicate real problems

**Slow processing**
- Reprojection is CPU-intensive for large images
- Downsize images before reprojection if needed: `image_small = image[::2, ::2]`

**Poor color balance**
- Try different scaling methods (`asinh`, `percentile`, etc.)
- Adjust the asinh strength: modify `a=0.02` in the code
- Check for bad pixels or cosmic rays before color combination

## Advanced: Custom Scaling

To modify scaling parameters:

```python
# Edit in color_image_notebook.py, line ~100
# Current: AsinhStretch(a=0.02)
# Try: AsinhStretch(a=0.05) for more aggressive stretch
```

Or subclass and override:

```python
class MyColorImageMaker(ColorImageMaker):
    def scale_image(self, image, method=None):
        # Custom scaling logic here
        return super().scale_image(image, method)
```
