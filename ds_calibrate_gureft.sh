#!/bin/bash
#
#SBATCH --job-name=DarkSage
#SBATCH --output=DS_calibrate_dump_12.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=bigmem
#SBATCH --time=1-00:00:00
#SBATCH --mem=500G

source /vast/palmer/apps/avx2/software/miniconda/24.7.1/etc/profile.d/conda.sh
conda activate DS24_env
module load GCC/12.2.0
module load GSL/2.5
module load OpenMPI/4.1.4-GCC-12.2.0
module load texlive/20220321-intel-2020b

cd /home/ajp79/project/fs1/DarkSage/priya_proj/runGUREFT90/DarkSage-astevens-disc/
python autocalibration/main.py -c input/gureft.par -b /home/ajp79/project/fs1/DarkSage/priya_proj/runGUREFT90/DarkSage-astevens-disc/darksage -o /home/ajp79/project/fs1/DarkSage/priya_proj/runGUREFT90/DarkSage-astevens-disc/autocalibration/01/ -x 'SMF_z6(7-11)*1,QLF_z6(43-48)*1' -S autocalibration/space_GUREFT.txt -s 16 -t 'chi2' -m 10

