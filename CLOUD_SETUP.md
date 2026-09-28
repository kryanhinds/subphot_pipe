# Claude Code cloud environment setup

How to run `subphot_pipe` offline in a Claude Code on the web session (Ubuntu
24.04 container), for example to rework the GRB 260310A photometry.

The cloud environment holds **no credentials**. `subphot_credentials_template.py`
leaves the Fritz token and LT-archive passwords blank on purpose, and the
environment's network policy blocks the astronomy services anyway. Every input
has to be local.

## What the SessionStart hook does

`.claude/hooks/session-start.sh` (registered in `.claude/settings.json`) only
runs in cloud sessions (`CLAUDE_CODE_REMOTE=true`). It:

1. `apt-get install`s SExtractor 2.28, SWarp 2.41 and PSFEx 3.24 and adds
   `sex` and `swarp` symlinks, because Ubuntu names the binaries
   `source-extractor` and `SWarp`.
2. Creates `.venv/` and installs `requirements.txt` into it. Ubuntu's patched
   system setuptools breaks the `sip-tpv` and `docopt` builds, so pip has to
   run in a venv.
3. Copies the credential-free template to `subphot_credentials.py` if that
   file is missing.
4. Writes `config_files/`. Every file comes from the pipeline's own
   generators except `sex.conv`, which is SExtractor's default 3x3 FWHM=2
   mask (the same filter `autoastrometry.py` writes).
5. Sets astropy's `auto_download = False` for IERS tables. They are used only
   for the sunset times printed at start-up.

## Local inputs to provide

| Put in the repo root | Used for |
|---|---|
| `data/<object>/*.fits` | science frames |
| `ref_imgs/...` (NOT references, passed with `-refimg`) | template images, so there is no PS1/Legacy/SDSS download |
| `ps_catalogs/ps_<ra>_<dec>_<rad>.xml` | cached PS1 catalogues; the whole folder from the machine that already ran the reduction |
| CFHT u-band catalogue, passed with `-refcat` | u-band zeropoint and catalogue fine-registration |

## Every outbound call in the pipeline and what triggers it

| Call (file) | Host | When it fires | Offline status |
|---|---|---|---|
| `Observer.at_site()` (`subphot_subtract.py`, `subphot_quicklook_pipe.py`) | astropy site list (`astropy.org` / `astropy.github.io`) | **every run**, at start-up. Sunset is used only to decide which night's date (`DATE`) names the output folders; photometry doesn't use it | **Fixed**: `observer_at_site()` falls back to built-in La Palma/Palomar coordinates, so the date logic is unchanged. Before this it raised `UnknownSiteException` offline and the run never started |
| astropy IERS tables | `datacenter.iers.org`, `maia.usno.navy.mil` | every run (sunset times) | disabled by the hook |
| `panstarrs_query` (`subphot_functions.py`) | `archive.stsci.edu` | g/r/i/z zeropoint, alignment and distortion catalogues, stack WCS check | cached: reads `ps_catalogs/ps_<ra>_<dec>_<rad>.xml` when present. On a cache miss it downloads, and offline the reduction crashes |
| `sdss_query` (`subphot_functions.py`) | `skyserver.sdss.org` | u band or `-s SDSS` / `-sdsscat`, **only when `-refcat` is `auto`** | not called when the CFHT u-band catalogue is passed with `-refcat`. Note it isn't cached |
| `make_sdss_ref` / `sdss_query_image` | `skyserver.sdss.org`, `dr16.sdss.org` | u band or SDSS references when `-refimg` is `auto` | avoided with `-refimg` |
| `panstamps` | `ps1images.stsci.edu` | g/r/i/z references when `-refimg` is `auto` and nothing is cached | avoided with `-refimg` |
| `download_legacy_survey_fits` | `www.legacysurvey.org` | `-s legacy` when `-refimg` is `auto`; skips the download if already in `ref_imgs/` | avoided with `-refimg` |
| `autoastrometry` | `tdc-www.harvard.edu`, `skyserver.sdss.org` | only with `-reastrom` on non-SEDM frames | don't pass `-reastrom` |
| Fritz `api` / `SN_data_phot` / photometry POST | `fritz.science` | only with `-up` / `-upf` | don't pass these; the token is blank |
| LT archive / quicklook downloads | `telescope.livjm.ac.uk` | only morning-roundup / `-qdl` download modes | not used; no passwords are present |

`astropy_ps1_astrometry.py` (MAST) is not imported by the v2 pipeline.
