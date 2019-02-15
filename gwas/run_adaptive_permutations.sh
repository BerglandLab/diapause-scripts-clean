#!/bin/bash

parameterFile="/scratch/$USER/genome-reconstruction/adaptive_permutation_input.txt"

pop=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $parameterFile | cut -f 1 )
pheno=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $parameterFile | cut -f 2 )
perm=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $parameterFile | cut -f 3 )
seed=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $parameterFile | cut -f 4 )

Rscript /scratch/pae3g/genome-reconstruction/scripts/genesis_adaptive_perm_3chr.R $pop $pheno $perm $seed
