#!/bin/bash
#
#SBATCH --job-name=meta-TCT-simulations
#SBATCH --ntasks=1 --cpus-per-task=2
<<<<<<< HEAD
#SBATCH --time=00:30:00
=======
#SBATCH --time=00:45:00
>>>>>>> a23c0d86bb3275fbca8da12859da6844315b8287
#SBATCH --cluster=wice
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -A lp_doctoralresearch

export OMP_NUM_THREADS=1

module use /apps/leuven/rocky8/icelake/2022b/modules/all
module load GSL
module load CMake
module load  R/4.2.2-foss-2022b

Rscript -e "renv::status()" -e "renv::restore()" - "renv::status()"




