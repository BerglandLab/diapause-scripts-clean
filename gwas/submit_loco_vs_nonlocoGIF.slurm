#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --mem 900G
#SBATCH -t 0-24:00:00
#SBATCH -p largemem
#SBATCH --account bergland-erickson
module load  gcc/7.1.0  openmpi/3.1.4

module load R/3.6.3

Rscript /scratch/pae3g/scripts/loco_vs_nonloco_GIF.R
