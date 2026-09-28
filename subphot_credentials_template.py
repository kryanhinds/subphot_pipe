"""Template for subphot_credentials.py (the real file is untracked).

Contains no credentials on purpose. Fritz, email and LT-archive access are
left blank, so -up / -e / morning-roundup downloads cannot authenticate from
an environment built from this template. Put real values only in a local,
untracked subphot_credentials.py on your own machines.
"""
import os
import shutil

_here = os.path.dirname(os.path.abspath(__file__)) + '/'

# Working directory for data, references and intermediate products (trailing slash)
path = os.environ.get('SUBPHOT_PATH', _here)
data1_path = os.environ.get('SUBPHOT_DATA1_PATH', path)

# External binaries (Ubuntu names SExtractor 'source-extractor' and SWarp 'SWarp')
sex_path = shutil.which('sex') or shutil.which('source-extractor') or 'sex'
swarp_path = shutil.which('swarp') or shutil.which('SWarp') or 'swarp'
psfex_path = shutil.which('psfex') or 'psfex'
panstamps_path = shutil.which('panstamps') or 'panstamps'
solve_field_path = shutil.which('solve-field') or 'solve-field'
solve_field_config_path = ''

# Deliberately blank: no Fritz token or LT archive passwords.
token = ''
proposals_arc = {}

# Pipeline tuning
image_size = 1500      # reference cutout size (pix)
starscale = 1.5        # detection threshold (bkg std)
search_rad = 1         # catalogue match radius (arcsec)
store_lc_ims = False
