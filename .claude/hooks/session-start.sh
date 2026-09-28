#!/bin/bash
# SessionStart hook for Claude Code on the web: installs the external
# astronomy binaries and Python dependencies subphot_pipe needs.
set -euo pipefail

if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

cd "$CLAUDE_PROJECT_DIR"

# SExtractor, SWarp and PSFEx from the Ubuntu archive
need_apt=0
for b in source-extractor swarp psfex; do command -v "$b" >/dev/null 2>&1 || need_apt=1; done
if [ "$need_apt" = 1 ]; then
  export DEBIAN_FRONTEND=noninteractive
  apt-get update -qq
  apt-get install -y -qq --no-install-recommends source-extractor swarp psfex >/dev/null
fi
# Pipeline code calls SExtractor as 'sex' in places; Ubuntu ships SWarp as 'SWarp'
if ! command -v sex >/dev/null 2>&1; then
  ln -sf "$(command -v source-extractor)" /usr/local/bin/sex
fi
command -v swarp >/dev/null 2>&1 || ln -sf "$(command -v SWarp)" /usr/local/bin/swarp

# Project virtualenv: Ubuntu's patched system setuptools breaks sdist builds
# (sip-tpv, docopt), so keep pip away from /usr/lib/python3/dist-packages
if [ ! -x .venv/bin/python ]; then
  python3 -m venv .venv
fi
.venv/bin/pip install -q --upgrade pip setuptools wheel
.venv/bin/pip install -q -r requirements.txt

# Per-machine config: env-var driven credentials (no secrets in the repo)
if [ ! -f subphot_credentials.py ]; then
  cp subphot_credentials_template.py subphot_credentials.py
fi
mkdir -p config_files temp_config_files ref_imgs data out

# config_files/ is untracked and host-specific. sex.conv is the one file the
# pipeline never writes itself: use SExtractor's default 3x3 FWHM=2 mask (the
# same filter autoastrometry.py writes). Generate the rest once here.
[ -f config_files/sex.conv ] || cp /usr/share/source-extractor/default.conv config_files/sex.conv
.venv/bin/python - <<'PY'
import os
from subphot_credentials import path
from subphot_functions import (writeswarpdefaultconfigfile, writeswarpconfigfile,
                               writesdssswarpconfigfile, writepsfexparfile,
                               prepsexfile, psfexfile)
import subphot_align_quick as saq
cf = lambda n: os.path.exists(path + 'config_files/' + n)
if not cf('config.swarp'): writeswarpdefaultconfigfile()
if not cf('config_comb.swarp'): writeswarpconfigfile()
if not cf('swarp_sdss.conf'): writesdssswarpconfigfile()
if not cf('default.param'): writepsfexparfile()
if not cf('prepsfex.sex'): prepsexfile()
if not cf('psfex_conf.psfex'): psfexfile()
if not cf('align_temp.param'): saq.writeparfile(path)
if not cf('align_sex.config'): saq.writeconfigfile()
PY

# No network in the sandbox: stop astropy trying to fetch IERS tables
# (maia.usno.navy.mil, datacenter.iers.org). Only sunset times use them; the
# bundled astropy-iers-data tables are ample for that.
mkdir -p "$HOME/.astropy/config"
if ! grep -qs '^\[utils.iers.iers\]' "$HOME/.astropy/config/astropy.cfg"; then
  cat >> "$HOME/.astropy/config/astropy.cfg" <<'CFG'

[utils.iers.iers]
auto_download = False
iers_degraded_accuracy = warn
CFG
fi

echo 'export PATH="'"$CLAUDE_PROJECT_DIR"'/.venv/bin:$PATH"' >> "$CLAUDE_ENV_FILE"
echo 'export PYTHONPATH="'"$CLAUDE_PROJECT_DIR"'${PYTHONPATH:+:$PYTHONPATH}"' >> "$CLAUDE_ENV_FILE"
