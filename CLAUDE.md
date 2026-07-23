# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`subphot_pipe` is an astronomical photometry pipeline for reducing multi-telescope optical imaging data. It performs image subtraction, aperture photometry, astrometry, and automated result upload to Fritz SkyPortal. Primary science targets are supernovae and other optical transients.

## Running the Pipeline

The main entry point is `subphot_subtract_v1.py`:

```bash
# Reduce all images in a folder (auto-download PS1 reference)
python subphot_subtract_v1.py -f data/ZTF25aceeneu/ -s PS1

# Reduce with SDSS references, upload to Fritz, and send email
python subphot_subtract_v1.py -f data/ZTF25aceeneu/ -s SDSS -up -e

# Reduce specific filters only
python subphot_subtract_v1.py -f data/ZTF25aceeneu/ -fb r i

# Reduce specific image files
python subphot_subtract_v1.py -i /path/to/image.fits

# Morning roundup (batch processing of last night's data with multiprocessing)
python subphot_subtract_v1.py -mrup -mp active

# Stack images where possible
python subphot_subtract_v1.py -f data/ZTF25aceeneu/ -stk

# Forced photometry at header RA/DEC
python subphot_subtract_v1.py -f data/ZTF25aceeneu/ -fp True

# Plot light curves
python subphot_subtract_v1.py -f data/ZTF25aceeneu/ -lc

# List FITS file metadata in a directory
python subphot_subtract_v1.py -ls data/ZTF25aceeneu/
```

Key flags:
- `-s` / `--survey`: Reference catalog — `PS1` (default Pan-STARRS) or `SDSS`
- `-tel` / `--telescope_facility`: `LT` (Liverpool Telescope, default), `HCT`, `SEDM`, etc.
- `-mp` / `--multipro`: Enable multiprocessing (`active`)
- `-cln` / `--cleandirs`: Clean intermediate products after reduction
- `-up` / `--upfritz`: Upload photometry to Fritz SkyPortal
- `-zp_only`: Compute and save zeropoint only (no subtraction)
- `-reastrom`: Force redo astrometry even if WCS is in header

## Configuration

All site-specific paths, API tokens, and credentials are in `subphot_credentials.py`:
- `path` / `data1_path`: Working directory for data and temporary files
- `token`: Fritz SkyPortal API token
- `swarp_path`, `sex_path`, `psfex_path`, `panstamps_path`: External binary paths
- `image_size`: Reference image cutout size in pixels (default 1500)
- `starscale`: Star detection threshold in units of background std (default 1.5)
- `search_rad`: Catalog matching radius in arcseconds (default 1)

## Architecture

### Core modules

- **`subphot_subtract_v1.py`** — CLI entry point; parses arguments, orchestrates single-image and batch reduction via `multi_subtract`
- **`subphot_quicklook_pipe_v1.py`** — Primary pipeline logic; contains `subtracted_phot` (single image) and `multi_subtract` (batch) classes; handles background subtraction, PSF measurement, template subtraction, photometry, zeropoint calibration, and Fritz upload
- **`subphot_functions.py`** — Utility functions: SWarp/SExtractor subprocess wrappers, PS1/SDSS catalog queries, astrometry helpers, email sending
- **`subphot_align_quick.py`** — Fast image alignment using SExtractor source extraction and star cross-matching
- **`subphot_align.py`** — Alternative alignment with star detection and coordinate transformation
- **`autoastrometry.py`** — Standalone astrometric plate solver (SExtractor-based, queries online catalogs)
- **`astropy_ps1_astrometry.py`** — Astrometric solver using Pan-STARRS catalog via astroquery
- **`psf_measure.py`** — PSF characterization via Gaussian fitting
- **`subphot_telescopes.py`** — `header_kw` dict mapping telescope names to FITS header keywords; add new telescopes here
- **`subphot_credentials.py`** — All configuration (paths, tokens, credentials)

### Processing flow

```
subphot_subtract_v1.py (CLI + arg parsing)
    └─ multi_subtract (subphot_quicklook_pipe_v1.py)
           └─ subtracted_phot (per image)
                  1. Read FITS, extract header metadata via subphot_telescopes.header_kw
                  2. Background estimation & subtraction (sigma-clipped)
                  3. Cosmic ray rejection (astroscrappy)
                  4. Astrometry: autoastrometry / astropy_ps1_astrometry
                  5. Download PS1/SDSS reference image (panstamps / astroquery)
                  6. Align science to reference (SWarp or subphot_align_quick)
                  7. PSF measurement on aligned image (psf_measure)
                  8. Image convolution to match PSFs
                  9. Template subtraction (reproject + direct subtraction)
                 10. Source extraction on subtracted image (SExtractor)
                 11. Aperture photometry at target position (photutils)
                 12. Zeropoint calibration via PS1/SDSS field stars
                 13. Save photometry CSV, optionally upload to Fritz
```

### External dependencies

These system binaries must be installed and paths set in `subphot_credentials.py`:
- **SExtractor** (`sex`) — source extraction
- **SWarp** — image resampling and co-addition
- **PSFEx** — PSF modelling (optional, `-psfex` flag)
- **panstamps** — PS1 image download

### Adding a new telescope

Add an entry to the `header_kw` dict in `subphot_telescopes.py` with the telescope name as key and a dict mapping standard field names (`filter`, `object`, `ra`, `dec`, `airmass`, `exptime`, `date`, `mjd`, `seeing`, `pixscale`, `gain`, etc.) to the corresponding FITS header keywords. Use `'-'` for unavailable fields and a numeric literal for fixed values (e.g. pixel scale).

### Data directories (runtime-generated)

- `data/<object_name>/` — science images and output photometry per target
- `data/zeropoints/by_obs_date/` — saved zeropoints
- `aligned_images/` — WCS-aligned science images
- `bkg_subtracted_science/` — background-subtracted images
- `convolved_psf/` — PSF-convolved images for subtraction
- `ref_imgs/` — cached reference images from PS1/SDSS
