# Overview

This set of scripts is used to generated a simulated F4 and F5 hybrid
population from the founders used in each of the two swarms. Simulated reads are then
produced at 0.5X coverage for each. These reads are mapped to produce
bam files, that then go through the reconstruction pipeline. In the
end, the actual sequence is compared to the reconstructed sequence to
determine the reconstruction accuracy.

## prepare files

* prepare\_read\_simulation\_fasta.sh takes founders vcf file, splits
  it into individual vcfs, then makes fasta files for each chromosome
  of each individual, using ambiguous bases for heterozygous
  sites.
* these files are used to produce the sequences of the simulated
  hybrid individuals

## generate simulated hybrid swarm

* nannyGenerateMapFreqs.slurm runs the GenerateMapFreqs.slurm script
  for different parameter combinations to generate the simulated
  swarms and produce the simulated mapped reads (bwa files) from forwardSimulator.Rscript
* There are 4 versions of this script for each combination of swarm A
  and B, generations 4 and 5
  * extract\_simulated_\_bams.sh pulls out the bam files that result from the simulation so they can be passed to the reconstruction
* summarize\_simulated\_haplotypes.R makes a combined file of all haplotypes produced by the simulation

## reconstruct simulated genomes

* extract\_simulated\_bams.sh extracts the bam files produced by the
simulation and renames them to include the simulation info in the filename

* reconstruct\_sim\_PAE.14.sh runs the reconstruction on simulated bams

* summarize\_simulated\_reconstruction.R reads in the reconstructed haplotypes into a single table

## calculate accuracy

* calc\_accuracy\_array.sh compares the genotypes of the raw
  reconstructed genomes to the known sequence from the simulated
  hybrid individual. The output is summarized with summarize\_accuracy.R

* accuracy is determined as the fraction of sites with identical
genotypes
