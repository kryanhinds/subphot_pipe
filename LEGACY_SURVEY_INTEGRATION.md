# Legacy Survey Integration

Legacy Survey (DECaLS) reference images are now integrated into the subphot_subtract.py pipeline, alongside PS1 and SDSS.

## Usage

### Basic Command

```bash
# Use Legacy Survey for g-band
cp3 subphot_subtract.py -f data/ZTF26aakjzdt_LT/20260420 -cut -o ZTF26aakjzdt -tel LT -fb g -s legacy

# Use Legacy Survey for multiple bands
cp3 subphot_subtract.py -f data/ZTF26aakjzdt_LT/20260420 -cut -o ZTF26aakjzdt -tel LT -fb g r -s legacy
```

### Command Options

```bash
# Download smaller cutout (default is 1200x1200 px)
-s legacy           # Use Legacy Survey for reference images

# Filter bands: g, r, i, z (DECaLS is optical only)
-fb g               # Single band
-fb g r i           # Multiple bands
```

## How It Works

1. When you specify `-s legacy`, the pipeline:
   - Detects the survey type
   - For each filter (g, r, i, z), downloads a 1200×1200 pixel cutout
   - Saves to `ref_imgs/legacysurvey_g.fits` (etc.)
   - Uses that as the reference image for alignment and photometry

2. If the download fails:
   - Pipeline automatically falls back to SDSS
   - Logs a warning message

## Advantages of Legacy Survey

✅ **Larger cutouts** — 1200×1200 px (vs your small 512×512)
✅ **No resizing cascade** — Avoids alignment problems
✅ **More calibration stars** — Better photometry constraints
✅ **No artifacts at your target** — Unlike PS1 which had an artifact at ZTF26aakjzdt

## Comparison of Reference Image Sources

| Feature | PS1 (panstamps) | SDSS | Legacy Survey |
|---------|---|---|---|
| Optical bands | g, r, i, z | g, r, i, z, u | g, r, i, z |
| Typical cutout size | ~300-500 px | Variable | 1200 px |
| Download method | panstamps CLI | astroquery | HTTP API |
| Coverage | ~3/4 sky | ~1/3 sky (declination < 27.5°) | ~3/4 sky |
| Availability | Fast | Fast | Fast |

## Examples

### Use case 1: Target has PS1 artifact

```bash
# Instead of PS1 (which has artifact)
cp3 subphot_subtract.py -f data/ZTF26aakjzdt_LT -s legacy -fb g
```

### Use case 2: Target outside SDSS footprint

```bash
# Your target is at Dec=+71.8°, outside SDSS southern limit
# Use Legacy Survey instead
cp3 subphot_subtract.py -f data/targets_north/ -s legacy -fb g r
```

### Use case 3: Default reference selection

The pipeline will choose references in this order:
1. Manual `-refimg` path (if provided)
2. Legacy Survey (if `-s legacy`)
3. PS1 (if `-s PS1` or default)
4. SDSS (if `-s SDSS` or u-band or SDSS-only mode)

## Troubleshooting

**Problem: "Legacy Survey download failed, falling back to SDSS"**
- Your coordinates are outside Legacy Survey footprint
- Check: https://www.legacysurvey.org/
- Fallback to SDSS: no action needed, pipeline handles it automatically

**Problem: Download too slow**
- First download can take ~30-60 seconds due to server processing
- Subsequent runs use cached files if in same directory
- Check your internet connection

**Problem: Image quality worse than expected**
- Verify image is 1200×1200 (not smaller due to alignment)
- Check alignment quality logs for NCC value (should be >0.3)
- Compare cutouts: reference quality at target location

## Integrating Custom Parameters

If you want to customize the download (size, pixscale, DR version):

Edit `subphot_quicklook_pipe.py` around line 2244:

```python
ref_path = download_legacy_survey_fits(
    ra=self.sci_c.ra.deg,
    dec=self.sci_c.dec.deg,
    band=self.sci_filt,
    output_dir=self.path+'ref_imgs/',
    size=1200,          # ← Change to 800, 1000, 1500, etc.
    pixscale=0.262,     # ← Keep at 0.262 for DECaLS native
    layer='ls-dr9',     # ← Change to 'ls-dr10' when available
    logger=self.sp_logger
)
```

## Files Modified

- `subphot_subtract.py` — Updated survey help text
- `subphot_quicklook_pipe.py` — Added Legacy Survey reference selection logic + survey attribute
- `subphot_functions.py` — Added `download_legacy_survey_fits()` function

## API Documentation

### `download_legacy_survey_fits()`

```python
from subphot_functions import download_legacy_survey_fits

ref_path = download_legacy_survey_fits(
    ra=219.317292,          # Right ascension (degrees)
    dec=71.841750,          # Declination (degrees)
    band='g',               # Filter: g, r, i, z
    output_dir='ref_imgs/', # Output directory
    size=1200,              # Cutout size in pixels
    pixscale=0.262,         # Pixel scale (arcsec/px)
    layer='ls-dr9',         # Data release version
    logger=None             # Optional logger object
)

# Returns: str (path to FITS file) or None (if failed)
```
