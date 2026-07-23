# Server / New-Machine Setup

This repo tracks **code only**. Two things are per-machine and deliberately
untracked — each host keeps its own copies and a `git pull` never touches them:

1. **`subphot_credentials.py`** — paths, API tokens, email credentials.
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

## One-time migration on minar (in-place branch switch)

The live installation at `/data/sedmdrp/sedmpy/subphot_pipe` switches to the
new `main` branch in place. Its previous state stays available on the
`server-live` branch as the rollback point. Because the directory doesn't
move, `subphot_credentials.py` and the localized `config_files/` keep working
unchanged, and the entry point keeps its name (`subphot_subtract.py`) so
cron/wrapper commands don't change.

```bash
cd /data/sedmdrp/sedmpy/subphot_pipe

# safety copies of the per-machine files, plus untracked files that would
# collide with files the new branch tracks
B=/data/sedmdrp/subphot_backup_$(date +%Y%m%d)
mkdir -p $B
cp -r config_files subphot_credentials.py $B/
mv subphot_quicklook_pipe_v1.py subphot_quicklook_pipe_old.py utils.ipynb build_psf.py $B/ 2>/dev/null

git fetch origin
git checkout main
# the switch removes files tracked on server-live but not on main —
# including the localized configs; restore them:
cp $B/config_files/*.sex $B/config_files/*.psfex $B/config_files/*.config \
   $B/config_files/*.conv $B/config_files/*.param $B/config_files/*.swarp \
   $B/config_files/*.conf config_files/
git config pull.ff only
```

If `git checkout main` refuses because other untracked files would be
overwritten, move the files it lists into `$B` and rerun it.

Rollback at any time: `git checkout server-live` and restore the configs from
`$B` the same way.

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
