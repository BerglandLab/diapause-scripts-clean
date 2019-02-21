#!/bin/bash

parameterFile="/scratch/$USER/genome-reconstruction/adaptive_permutation_input2.txt"

pop=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $parameterFile | cut -f 1 )
pheno=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $parameterFile | cut -f 2 )
perm=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $parameterFile | cut -f 3 )
seed=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $parameterFile | cut -f 4 )

module load gcc
module load R/3.5.1

Rscript /scratch/pae3g/genome-reconstruction/scripts/genesis_adaptive_perm_3chr.R $pop $pheno $perm $seed

#Use this to run 500 permutations of selected snps
#sbatch --array=1-1%100 --account=bergland-erickson --ntasks-per-node=20 --partition=standard --mem=192G --time=0-10:00:00 /scratch/pae3g/genome-reconstruction/scripts/run_adaptive_permutations2.sh
