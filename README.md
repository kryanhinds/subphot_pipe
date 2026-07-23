# SEDM Photometry Pipeline

Astronomical photometry pipeline for reducing multi-telescope optical imaging data.

## Quick Start

```bash
# Reduce all images in a folder
python subphot_subtract_v1.py -f data/ZTF25abnjzp/ -s PS1

# With upload to Fritz and email
python subphot_subtract_v1.py -f data/ZTF25abnjzp/ -s PS1 -up -e

# Morning roundup with multiprocessing
python subphot_subtract_v1.py -mrup -mp active
```

See `CLAUDE.md` for full documentation.
