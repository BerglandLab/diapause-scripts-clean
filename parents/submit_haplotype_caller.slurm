#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem 64000
#SBATCH -t 0-24:00:00
#SBATCH -p standard
#SBATCH --account bergland-erickson

/scratch/pae3g/parents/scripts/haplotype_caller_rivanna.sh "$1"
