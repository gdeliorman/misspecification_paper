#!/bin/bash
#
#SBATCH --job-name=meta-TCT-simulations
#SBATCH --ntasks=1 --cpus-per-task=72
#SBATCH --time=20:00:00
#SBATCH --cluster=wice
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -A lp_doctoralresearch

export OMP_NUM_THREADS=1

module use /apps/leuven/rocky8/icelake/2022b/modules/all
module load GSL
module load CMake
module load  R/4.2.2-foss-2022b

Rscript -e "renv::restore()" -e "Sys.setenv(TZ='Europe/Brussels')"
Rscript real-data-sensitivity-analysis.R 71



