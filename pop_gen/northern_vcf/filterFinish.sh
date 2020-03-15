#! /bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=64000
#SBATCH --time=48:00:00
#SBATCH --partition=economy
#SBATCH --account=berglandlab

### filter sites based on PASS field
/home/ab5dr/vcftools/bin/vcftools \
--vcf /scratch/ab5dr/InbredLines/filter/BME_LN_LNPA.raw.vsqr.vcf \
--remove-filtered-all \
--recode \
--recode-INFO-all \
--out /scratch/ab5dr/InbredLines/filter/BME_LN_LNPA.raw.vsqr.filter
	
### remove repetative regions & filter around indels & convert to BCF
/scratch/ab5dr/bedtools/bedtools2/bin/bedtools intersect -sorted -v -header \
-a /scratch/ab5dr/InbredLines/filter/BME_LN_LNPA.raw.vsqr.filter.recode.vcf \
-b /scratch/ab5dr/InbredLines/filter/repMasker.dgrp2_indel_100bp.sort.bed | \
tee /scratch/ab5dr/InbredLines/filter/BME_LN_LNPA.raw.vsqr.filter.recode.noREP.noINDEL.vcf | \
/home/ab5dr/htslib/bin/bgzip -c > \
/scratch/ab5dr/InbredLines/filter/BME_LN_LNPA.raw.vsqr.filter.recode.noREP.noINDEL.vcf.gz

### tabix
/home/ab5dr/htslib/bin/tabix -p vcf /scratch/ab5dr/InbredLines/filter/BME_LN_LNPA.raw.vsqr.filter.recode.noREP.noINDEL.vcf.gz
