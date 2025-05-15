import argparse
import sys
from subphot_credentials import *
import os

parser = argparse.ArgumentParser()
parser.add_argument('--folder','-f',default='', nargs='*',
    help='List of folders to reduce (accepts wildcards or '
    'space-delimited list)')

parser.add_argument('--multipro','-mp',default='inactive',
    help="Speed up subtraction by utilising multiprocessing and using 5 nodes")

parser.add_argument('--bands','-fb',default=['All'],nargs='+',
    help="Filters to process, default is All e.g. g,r,i,z,u for SDSS-G, SDSS-R, SDSS-I, SDSS-Z & SDSS-U")

parser.add_argument('--sci_names','-sn',default=['All'],nargs='+',
    help="Science Names - names of science object to process, default is all found")

parser.add_argument('--out_dir','-o',default='by_name',
    help="Dest. for output photometry files, default creates a folder for photometry by \
    observation date e.g. for h_e_20220206_***_1_1.fits output photometry dest will be photometry_data/20220206")

parser.add_argument('--upload','-up',action='store_true',default=False,
    help="Upload photometry to Fritz (SkyPortal), default is False. If True, will require correct credentials in subphot_credentials.py (Fritz token enabled for uploading)")

parser.add_argument('--upf','-upf',action='store_true',default=False,
    help="Force upload photometry to Fritz (SkyPortal), default is False. If True, will require correct credentials in subphot_credentials.py (Fritz token enabled for uploading)")

parser.add_argument('--sedm','-sedm',default=False,action='store_true',help="Telescope is SEDM")

args = parser.parse_args()

if not os.path.exists(path+'subphot_shell_scripts'):os.mkdir(path+'subphot_shell_scripts')
if args.sedm==True:sedm_tel='-tel sedm'
else: sedm_tel=''

if args.upload==True:upload = '-up'
else: upload = ''

if args.sci_names != 'All':sn = '-sn '+' '.join(args.sci_names)
else: sn = ''

if args.bands[0] != 'All':fb = '-fb '+' '.join(args.bands)
else: fb = ''

if args.upf==True:upf = '-upf'
else: upf = ''

print('Submitting jobs for folders: ',' '.join(args.folder))
for fold in args.folder:
    if '/' in fold:fold_log = fold.split('/')[-1]
    else: fold_log = fold
    fold_file_ = f'''#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --threads-per-core=1
#SBATCH --partition=transients
#SBATCH --job-name={fold}
#SBATCH --output=/mnt/data1/users/arikhind/phot_data/phot_logs/{fold_log}.log


flight env activate gridware
bash
module load apps/swarp/2.41.5/gcc-8.5.0 apps/psfex/3.24.1/gcc-8.5.0 apps/sextractor/2.28.0/gcc-8.5.0
conda activate subphot


python3 ~/subphot_pipe/subphot_subtract.py -f {fold} -o {args.out_dir} -swarp -psfex -log {fold}_log_t -cut {sedm_tel} {upload} {sn} {fb} {upf}
    '''

    if '/' in fold:fold = fold.split('/')[-1]
    with open(path+f'subphot_shell_scripts/subphot_rf_{fold}.sh','w') as f:
        f.write(fold_file_)
        f.close()

    print('Written ',path+f'subphot_shell_scripts/subphot_rf_{fold}.sh')

submit = input('Press Enter to submit jobs')
if submit == '':
    for fold in args.folder:
        if '/' in fold:fold = fold.split('/')[-1]
        os.system(f'sbatch {path}subphot_shell_scripts/subphot_rf_{fold}.sh')
        print(f'Submitted {fold}')
        

