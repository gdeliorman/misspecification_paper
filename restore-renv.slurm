#!/bin/bash
#
#SBATCH --job-name=restore-renv-misspecification
#SBATCH --ntasks=1 --cpus-per-task=2
#SBATCH --time=01:45:00
#SBATCH --cluster=wice
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -A lp_doctoralresearch

export OMP_NUM_THREADS=1

module use /apps/leuven/rocky8/icelake/2022b/modules/all
module load GSL
module load CMake
module load  R/4.2.2-foss-2022b

Rscript -e "renv::status()" -e "renv::restore()" - "renv::status()"




