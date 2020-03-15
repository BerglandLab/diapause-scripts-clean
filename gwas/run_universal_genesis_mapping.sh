#!/bin/bash

parameterFile="/scratch/$USER/genome-reconstruction/universal_input.txt"


pop=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $parameterFile | cut -f 1 )
phenotype=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $parameterFile | cut -f 2 )
draw=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $parameterFile | cut -f 3 )
perm=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $parameterFile | cut -f 4 )
seed=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $parameterFile | cut -f 5 )
model=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $parameterFile | cut -f 6 )

module load gcc
module load R/3.5.1

Rscript /scratch/pae3g/scripts/universal_genesis_mapping.R $pop $phenotype $draw $perm $seed $model

# sbatch --array=1-10%100 --account=bergland-erickson --ntasks-per-node=20 --partition=standard --mem=150G --time=0-4:00:00 /scratch/$USER/scripts/run_universal_genesis_mapping.sh
