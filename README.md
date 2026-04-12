# SEDM Photometry Pipeline

Astronomical photometry pipeline for reducing multi-telescope optical imaging data.

## Quick Start

```bash
# Reduce all images in a folder
python sedm_subtract.py -f data/ZTF25abnjzp/ -s PS1

# With upload to Fritz and email
python sedm_subtract.py -f data/ZTF25abnjzp/ -s PS1 -up -e

# Morning roundup with multiprocessing
python sedm_subtract.py -mrup -mp active
```

See `CLAUDE.md` for full documentation.
