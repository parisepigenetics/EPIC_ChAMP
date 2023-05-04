#!/bin/bash
################################ Slurm options #################################
### Job name
#SBATCH --job-name=EPIC
### Output
#SBATCH --output=EPIC-%j.out  # both STDOUT and STDERR
##SBATCH -o slurm.%N.%j.out  # STDOUT file with the Node name and the Job ID
##SBATCH -e slurm.%N.%j.err  # STDERR file with the Node name and the Job ID
### Limit run time "days-hours:minutes:seconds"
##SBATCH --time=24:00:00
### Requirements
#SBATCH --partition=ipop-up
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=15GB
### Email
##SBATCH --mail-user=email@address
##SBATCH --mail-type=ALL
################################################################################



module purge
module load snakemake

srun snakemake --use-conda --conda-frontend conda --core 5 --jobs=1 -s ChAMP.rules

mkdir output 
mv *.out output
rm pD_ChAMP.csv

