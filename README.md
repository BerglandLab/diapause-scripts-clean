# Introduction

This project contains all scripts necessary for the analysis of data for Erickson et al (submitted), examining the genetic basis of diapause in the model fruit fly _Drosophila melanogaster_.

All sequence data has been deposited in the SRA in BioProject # PRJNA522357

Most data files necessary for analysis and plotting are available on DataDryad.

The general workflow procedes through the folders in the following order:

1. **parents** contains all scripts for the sequencing and variant analysis of the 68 parental strains. These scripts produce the final vcf files used for reconstruction of hybrid genotypes

2. **hybrid_reconstruction** contains all scripts for processing the raw sequencing reads of hybrid individuals, reconstructing full diploid genotypes, and cleaning up these genotypes into a set of 100 vcfs with missing data imputed

3. **hybrid_simulation** contains scripts for generating simulated hybrid swarm genomes and reconstructing these genomes to estimate the accuracy of reconstruction with our experimental design

4.  **ld** contains scripts to analyze linkage disequilibrium in the hybrid swarm relative to parental populations and the DGRP

5. **gwas** contains the scripts for performing genome-wide association mapping of the diapause phenotype using the GENESIS software package in R. The scripts include the permutation analysis.

6. **clinal\_seasonal** contains the scripts used for looking at the enrichment of the diapause GWAS results relative to known clinal and seasonal SNPs from the data in Bergland and Machado et al 2019.

7. **pop\_gen** contains scripts for the  African allele frequency, and IHS analysis of GWAS SNPs

8. **other** contains several miscellaneous analysis items related to the phenotyping and field experiments

The final figures and some of the analysis are done and produced in the file "code\_by\_figure_R2.R".
