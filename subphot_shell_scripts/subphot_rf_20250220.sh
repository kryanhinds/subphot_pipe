flight env activate gridware
bash
module load apps/swarp/2.41.5/gcc-8.5.0 
module load apps/psfex/3.24.1/gcc-8.5.0 
module load apps/sextractor/2.28.0/gcc-8.5.0
conda activate subphot


python3 /users/arikhind/subphot_pipe/subphot_subtract.py -f Quicklook/20250220 | tee -a /users/arikhind/subphot_pipe/scront_test_20250220.log
