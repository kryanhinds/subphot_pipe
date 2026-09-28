"""Template for subphot_credentials.py (the real file is untracked).

Copy to subphot_credentials.py and edit, or leave it as-is and supply values
through environment variables. Never commit tokens or passwords: in the Claude
Code cloud environment set them as environment variables in the environment
settings instead (see CLOUD_SETUP.md).
"""
import os
import shutil

_here = os.path.dirname(os.path.abspath(__file__)) + '/'

# Working directory for data, references and intermediate products (trailing slash)
path = os.environ.get('SUBPHOT_PATH', _here)
data1_path = os.environ.get('SUBPHOT_DATA1_PATH', path)

# External binaries (Ubuntu names source-extractor 'source-extractor', not 'sex')
sex_path = os.environ.get('SEX_PATH') or shutil.which('sex') or shutil.which('source-extractor') or 'sex'
swarp_path = os.environ.get('SWARP_PATH') or shutil.which('swarp') or shutil.which('SWarp') or 'swarp'
psfex_path = os.environ.get('PSFEX_PATH') or shutil.which('psfex') or 'psfex'
panstamps_path = os.environ.get('PANSTAMPS_PATH') or shutil.which('panstamps') or 'panstamps'
solve_field_path = os.environ.get('SOLVE_FIELD_PATH') or shutil.which('solve-field') or 'solve-field'
solve_field_config_path = os.environ.get('SOLVE_FIELD_CONFIG_PATH', '')

# Fritz SkyPortal API token (only needed for -up / Fritz queries)
token = os.environ.get('FRITZ_TOKEN', '')

# Liverpool Telescope archive proposals: {proposal: [password, archive index, cat index]}
# Only needed for LT morning-roundup downloads; leave empty otherwise.
proposals_arc = {}

# Email (-e flag)
email_user = os.environ.get('SUBPHOT_EMAIL_USER', '')
email_password = os.environ.get('SUBPHOT_EMAIL_PASSWORD', '')

# Pipeline tuning
image_size = int(os.environ.get('SUBPHOT_IMAGE_SIZE', 1500))   # reference cutout size (pix)
starscale = float(os.environ.get('SUBPHOT_STARSCALE', 1.5))     # detection threshold (bkg std)
search_rad = float(os.environ.get('SUBPHOT_SEARCH_RAD', 1))     # catalogue match radius (arcsec)
store_lc_ims = os.environ.get('SUBPHOT_STORE_LC_IMS', '0') == '1'
