#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH --mem 100G
#SBATCH -t 0-72:00:00
#SBATCH -p standard
#SBATCH --account bergland-erickson
module load  gcc/7.1.0  openmpi/3.1.4

module load R/3.5.3


Rscript /scratch/pae3g/scripts/ld_decay_dgrp_maf0.05.R
