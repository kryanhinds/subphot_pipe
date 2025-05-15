#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --threads-per-core=1
#SBATCH --partition=transients
#SBATCH --job-name=ZTF24abbenwl
#SBATCH --output=/mnt/data1/users/arikhind/phot_data/phot_logs/ZTF24abbenwl.log


flight env activate gridware
bash
module load apps/swarp/2.41.5/gcc-8.5.0 apps/psfex/3.24.1/gcc-8.5.0 apps/sextractor/2.28.0/gcc-8.5.0
conda activate subphot


python3 ~/subphot_pipe/subphot_subtract.py -f ZTF24abjjpbo -o run_unsubtract_ZTF24abbenwl -swarp -psfex -log run_unsubtract_ZTF24abbenwl -cut -un
    