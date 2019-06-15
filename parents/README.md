# Overview

These scripts process new sequencing data for the inbred lines used
to initiate the hybrid swarms. The new sequencing data was combined
with existing sequencing data (available on the SRA) for variant
calling. Two vcfs are produced: one containing the most
high-confidence variants to be used for reconstruction, and a second
vcf with more variants used for association mapping.

## Download and process existing parental sequence reads

* Use sraData.delim and download\_sraData.sh to download appropriate
sequence files.
* Use map\_sra2.sh to map and clean up the
reads.

##  map parental reads sequenced here

* pear\_bwa.sh

## merge existing data and new data bam files


## genotype parents

* haplotype\_caller\_workflow.sh has all commands for executing calling haplotypes and making VCFs in GATK. This script uses haplotype\_caller\_rivanna.sh and submit\_haplotype\_caller.slurm
* filter\_parents.R does additional filtering of SNPs on missing parental snp information
*
