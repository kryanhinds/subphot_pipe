#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --threads-per-core=1
#SBATCH --partition=transients
#SBATCH --job-name=ZTF22aapubuy
#SBATCH --output=ZTF22aapubuy.log


flight env activate gridware
bash
module load apps/swarp/2.41.5/gcc-8.5.0 apps/psfex/3.24.1/gcc-8.5.0 apps/sextractor/2.28.0/gcc-8.5.0
conda activate subphot


python3 ~/subphot_pipe/subphot_subtract.py -f ZTF22aapubuy -o sedm_comps2 -swarp -psfex -log ZTF22aapubuy_log_t -cut -tel sedm
# python3 ~/subphot_pipe/subphot_subtract.py -f ZTF24aapvieu -o sedm_comps2 -swarp -psfex -log ZTF24aapvieu_log_t -cut -tel sedm



