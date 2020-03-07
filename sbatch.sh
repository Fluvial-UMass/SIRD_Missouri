#!/bin/bash
#SBATCH -n 20  # Number of Cores per Node
#SBATCH --mem=8192  # Requested Memory
#SBATCH -p cpu1  # Partition
#SBATCH -t 20:00:00  # Job time limit
#SBATCH -o dassim.out  # %j = job ID

conda activate geo
source /home/yishitsuka_umass_edu/.bashrc
python /project/pi_cjgleason_umass_edu/yuta/ISRD/MSR/src/dassim.py -c /project/pi_cjgleason_umass_edu/yuta/ISRD/MSR/src/config/config_MERIT_cal10.ini -s 20021001 -e 20101231