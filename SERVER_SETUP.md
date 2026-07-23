# Server / New-Machine Setup

This repo tracks **code only**. Two things are per-machine and deliberately
untracked — each host keeps its own copies and a `git pull` never touches them:

1. **`sedm_credentials.py`** — paths, API tokens, email credentials.
   The `path` variable must point at the pipeline directory (trailing slash)
   and the external binary paths (`sex_path`, `swarp_path`, `psfex_path`,
   `panstamps_path`) must match the host.
2. **`config_files/`** — SExtractor / SWarp / PSFEx configs. These embed
   absolute paths (`PARAMETERS_NAME`, `FILTER_NAME`, `PSF_DIR`,
   `CHECKIMAGE_NAME`, ...) and must be localized once per machine.
   The top-level `sex.config` / `sex.conv` are generated automatically by
   `autoastrometry.py` and need no setup.

Required files in `config_files/` (copy from an existing installation, then
fix the embedded paths):

```
align_sex.config  align_sex.conv  align_temp.param
config.swarp      config_comb.swarp  config_resize.swarp  swarp_sdss.conf
default.param     temp.param
prepsfex.sex      psfex_conf.psfex
sex.config        sex.conv
```

## One-time migration on minar (from the old subphot_pipe checkout)

The old installation lives at `/data/sedmdrp/sedmpy/subphot_pipe` (branch
`server-live` = snapshot of the live March 2026 state). Leave it untouched as
a fallback. Set up the new clone beside it, not inside `sedmpy`:

```bash
cd /data/sedmdrp
git clone https://github.com/kryanhinds/sedm-phot.git sedm_phot
cd sedm_phot
git config pull.ff only

# per-machine files from the old installation
OLD=/data/sedmdrp/sedmpy/subphot_pipe
cp $OLD/subphot_credentials.py sedm_credentials.py
mkdir -p config_files
cp $OLD/config_files/*.sex $OLD/config_files/*.psfex $OLD/config_files/*.config \
   $OLD/config_files/*.conv $OLD/config_files/*.param $OLD/config_files/*.swarp \
   $OLD/config_files/*.conf config_files/

# localize: configs and credentials still point at the old directory
sed -i "s|$OLD|/data/sedmdrp/sedm_phot|g" config_files/* sedm_credentials.py
```

Then edit `sedm_credentials.py` by hand and verify `path` / `data1_path` point
at `/data/sedmdrp/sedm_phot/`, and repoint whatever launches the pipeline
(cron entry / wrapper script) from the old directory to the new one.

## Routine updates (the whole point)

Develop and test on the local machine, then:

```bash
git push            # local
git pull            # on minar — fast-forward only; fails loudly if the
                    # server tree was edited, which should never happen
```

Never edit tracked files on the server. If a server-side hot-fix is
unavoidable, commit it to a branch and push it so it can be merged properly:

```bash
git checkout -b hotfix-YYYYMMDD && git commit -am "describe fix" && git push -u origin hotfix-YYYYMMDD
```
