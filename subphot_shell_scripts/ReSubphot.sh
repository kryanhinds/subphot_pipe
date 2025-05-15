#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --threads-per-core=1
#SBATCH --partition=transients
#SBATCH --job-name=RecentReanalysis
#SBATCH --output=/mnt/data1/users/arikhind/phot_data/phot_logs/RecentReanalysis.log


flight env activate gridware
bash
module load apps/swarp/2.41.5/gcc-8.5.0 apps/psfex/3.24.1/gcc-8.5.0 apps/sextractor/2.28.0/gcc-8.5.0
conda activate subphot


# python3 ~/subphot_pipe/subphot_subtract.py -f ZTF24abjjpbo -o run_unsubtract_ZTF24abjjpbo -swarp -psfex -log run_unsubtract_ZTF24abjjpbo -cut -un

python3 ~/subphot_pipe/subphot_subtract.py -f entire_lc_imgs/SN2024abfl -swarp -psfex -up

python3 ~/subphot_pipe/subphot_subtract.py -f entire_lc_imgs/SN2024abup -swarp -psfex -up

python3 ~/subphot_pipe/subphot_subtract.py -f Quicklook/20241124 -swarp -psfex -sn ZTF24abdiwwv -upf