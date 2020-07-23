#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=20
#SBATCH --mem 200G
#SBATCH -t 0-48:00:00
#SBATCH -p standard
#SBATCH --account bergland-erickson
module load  gcc/7.1.0  openmpi/3.1.4

module load R/3.5.3

Rscript /scratch/pae3g/scripts/sign_test_universal_threshold_dropmissing_droptlk.R
