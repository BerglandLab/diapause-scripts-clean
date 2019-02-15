#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem 6000
#SBATCH -t 0-48:00:00
#SBATCH -p standard
#SBATCH --account bergland-erickson

echo $HOSTNAME
echo $$

# Nanny script 

topDirectory="/scratch/$USER/genome-reconstruction/"
# change reference group to whatever directory this will reference; cageA, cageB, dgrp2, etc
referenceGroup=cageB
referenceGenome="$topDirectory"/variants/$referenceGroup/ref/all_dmel.fasta
nGenerations="5"
nFounders="34"
populationSize="10000"
# growthRate of 1 is constant population size
growthRate="1"
bedfile=/"$topDirectory"/variants/$referenceGroup/recombination.bed
# sex can be hermaphroditic or dioecious
sex="dioecious"
# recombination can be 'both' or 'femaleOnly'
recombination="femaleOnly"
# number of threads to use; for efficient use of cpu hours, should be a factor of the number of individuals done per population, and probably no greater than 5
nThreads="5"
# point to a vcf file, which is used to sample potential founders (only extracts #CHROM row with founder IDs)
vcf="$topDirectory"/variants/$referenceGroup/vcf/$referenceGroup.2L.vcf
# number of individuals to map per forward-simulated population
nIndsPerPop="10"
readLength=150
coverage=0.5
maxNjobs=4
nIndependentPopulations=100

cd $topDirectory

# infinite loop that checks check number of currently running jobs
for rep in $(seq 1 $nIndependentPopulations); do
    while true; do
        # count jobs running or pending
        nRUNNINGjobs=$(sacct -u $USER | grep 'RUNNING\|PENDING' | wc -l)

        # submit job IF below threshold
        if [[ "$nRUNNINGjobs" -lt "$maxNjobs" ]]; then
            sbatch /scratch/$USER/genome-reconstruction/scripts/GenerateMapFreqs.slurm \
                $topDirectory \
                $referenceGroup \
                $referenceGenome \
                $nGenerations \
                $nFounders \
                $populationSize \
                $growthRate \
                $bedfile \
                $sex \
                $recombination \
                $nThreads \
                $vcf \
                $nIndsPerPop \
                $readLength \
                $coverage \
                $rep
                
            
            # wait before submitting next job
            sleep 20s
            break
        fi
       
        # wait before checking running job count again
        sleep 5s
    done
done
