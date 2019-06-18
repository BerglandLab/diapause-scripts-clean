# Introduction

This project contains all scripts necessary for the analysis of data for Erickson et al (submitted), examining the genetic basis of diapause in the model fruit fly _Drosophila melanogaster_.

All sequence data has been deposited in the SRA in BioProject # PRJNA522357

Most data files necessary for analysis and plotting are available on DataDryad.

The workflow procedes through the folders in the following order:

1. **parents** contains all scripts for the sequencing and variant analysis of the 68 parental strains. These scripts produce the final vcf files used for reconstruction of hybrid genotypes

2. **hybrid_reconstruction** contains all scripts for processing the raw sequencing reads of hybrid individuals, reconstructing full diploid genotypes, and cleaning up these genotypes into a set of 100 vcfs with missing data imputed

3. **hybrid_simulation** contains scripts for generating simulated hybrid swarm genomes and reconstructing these genomes to estimate the accuracy of reconstruction with our experimental design

4. **gwas** contains the scripts for performing genome-wide association mapping of the diapause phenotype using the GENESIS software package in R. The scripts include the permutation analysis and the calculation of empirical p-values relative to permutation.

5. **clinal\_seasonal\_enrichment** contains the scripts used for looking at the enrichment of the diapause GWAS results relative to known clinal and seasonal SNPs from the data in Bergland and Machado et al 2019.

6. **pop\_gen** contains scripts for the Tajima's D, African allele frequency, and IHS analysis of index SNPs


The final figures are produced in the file "code\_by\_figure.R".
