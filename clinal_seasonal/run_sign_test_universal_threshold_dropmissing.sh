#!/bin/bash
module load  gcc/7.1.0  openmpi/3.1.4

module load R/3.5.3


Rscript /scratch/pae3g/scripts/sign_test_universal_threshold_dropmissing.R

# sbatch --account=bergland-erickson --ntasks-per-node=16 --partition=largemem --mem=200G --time=0-6:00:00 /scratch/$USER/scripts/run_sign_test_universal_threshold_dropmissing.sh
