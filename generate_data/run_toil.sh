#!/bin/sh
#SBATCH --job-name=molmimic
#SBATCH --mem=12000
#SBATCH --time=7-00:00:00
#SBATCH --partition=standard
#SBATCH -A muragroup
export TOIL_SLURM_ARGS="-t 10:00:00 --partition=standard -A muragroup"
python /project/ppi_workspace/molmimic/generate_data/build_full_dataset.py --defaultCores=4 --defaultMemory=12000M --maxLocalJobs=10000 --batchSystem=slurm default 
