# Claude Code cloud environment setup

How to run `subphot_pipe` in a Claude Code on the web session (Ubuntu 24.04
container), for example to rework the GRB 260310A / ZTF26aakjzdt SEDM
photometry.

## What the SessionStart hook does

`.claude/hooks/session-start.sh` (registered in `.claude/settings.json`) runs
at the start of every cloud session. It does nothing on local machines, where
`CLAUDE_CODE_REMOTE` is unset. It:

1. `apt-get install`s SExtractor 2.28, SWarp 2.41 and PSFEx 3.24 from the
   Ubuntu archive, and adds `sex` and `swarp` symlinks because Ubuntu names the
   binaries `source-extractor` and `SWarp`.
2. Creates `.venv/` and installs `requirements.txt` into it. Ubuntu's patched
   system setuptools breaks the sdist builds of `sip-tpv` and `docopt`, so pip
   has to run in a venv. `.venv/bin` is put first on `PATH` for the session.
3. Copies `subphot_credentials_template.py` to `subphot_credentials.py` if
   that file is missing. The template reads everything from environment
   variables and finds the binaries on `PATH`.
4. Writes `config_files/`. Every file comes from the pipeline's own
   generators except `sex.conv`, which is SExtractor's default 3x3 FWHM=2
   mask (the same filter `autoastrometry.py` writes).

The hook is idempotent. It takes about 75 s on a fresh container and about
13 s on a warm one.

## Things you must set in the environment settings

Open the cloud environment menu in the session title bar and choose **Edit**.

### 1. Network access (required)

With the default network policy, the proxy refuses every astronomy service
the pipeline uses. Add these hosts to the allowed domains, or choose a broader
access level:

| Host | Used for |
|---|---|
| `www.legacysurvey.org` | Legacy Survey reference cutouts (`-s legacy`, which the GRB run used) |
| `ps1images.stsci.edu` | PS1 reference images |
| `catalogs.mast.stsci.edu`, `archive.stsci.edu`, `mast.stsci.edu` | PS1 catalogue for zeropoints and astrometry |
| `skyserver.sdss.org`, `dr16.sdss.org`, `dr12.sdss.org` | SDSS catalogue and images (`-s SDSS`) |
| `tdc-www.harvard.edu` | `autoastrometry` catalogue queries |
| `fritz.science` | Fritz upload and queries (`-up`) |
| `datacenter.iers.org`, `hpiers.obspm.fr` | astropy IERS table updates (optional) |

PyPI, `archive.ubuntu.com` and conda-forge are already reachable.

### 2. Secrets (optional)

Set these as environment variables. Never commit them.

- `FRITZ_TOKEN`: the Fritz SkyPortal API token. Only needed for `-up`.
- `SUBPHOT_EMAIL_USER` and `SUBPHOT_EMAIL_PASSWORD`: only needed for `-e`.

Tuning overrides: `SUBPHOT_IMAGE_SIZE`, `SUBPHOT_STARSCALE`,
`SUBPHOT_SEARCH_RAD`, `SUBPHOT_PATH`.

## Getting the GRB data into the session

The science frames are not in git (`data/` and `grb*/` are ignored). To make
them available, pick one of these:

- Put `data/ZTF26aakjzdt/` (the SEDM `rc*_ZTF26aakjzdt_*.fits` frames) and,
  optionally, `ref_imgs/ZTF26aakjzdt_legacysurvey_*.fits` in Google Drive and
  ask Claude to fetch them through the Drive connector.
- Push them to a separate data branch or release asset.

Then run:

```bash
python subphot_subtract.py -f data/ZTF26aakjzdt/ -s legacy -tel SEDM
```
