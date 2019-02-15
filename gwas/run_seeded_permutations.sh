#!/bin/bash

parameterFile="/scratch/$USER/genome-reconstruction/seeded_permutation_input.txt"

file=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $parameterFile | cut -f 1 )
pop=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $parameterFile | cut -f 2 )
phenotype=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $parameterFile | cut -f 3 )
draw=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $parameterFile | cut -f 4 )
perm=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $parameterFile | cut -f 5 )
seed=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $parameterFile | cut -f 6 )


Rscript /scratch/pae3g/genome-reconstruction/scripts/genesis_loco_byswarm_final2_draw_seed.R $file $pop $phenotype $draw $perm $seed
