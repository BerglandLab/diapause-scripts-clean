# Overview

This set of scripts has everything required to perform association
mapping using the R package "GENESIS". GENESIS requires a genetic
relatedness matrix as a covariate. The association mapping is
performed on each of the 100 randomly drawn genotype VCF files, and
100 permutations of the phenotypes are performed. For each
permutation, the p-values are averaged across all 100
imputations. Empirical p-values are calculated by grouping SNPs into
bins based on allele frequency and missing data. Actual p values are
compared to the distribution of all permuted p-values for SNPs in the
same category.

## generate master info about each SNP

* hwe\_replace\_filters.R" calculates a number of statistics for each
SNP and uses these to generate filters for which SNPs will go into the
association mapping, and which SNPs need additional permutations

## make GDS and compute GRM for each imputation

* use "makeGDS\_LOCO\_GRM.R"

## run gwas for 100 imputations of raw data + 100 permutations

* generate input (file names, permutation numbers, and random seeds)
with "seeded\_permutation\_input.R"

 * use "run\_seeded\_permutations.sh" to run "genesis\_mapping.R" for
 all imputations/permutations

## summarize imputations

* use "summarize\_impuations.R" to calculate mean and median gwas statistics for
  all permutations

## perform adaptive permutations

* input is also generated in "seeded\_permutation\_input.R"
* use "run\_adaptive\_permutations.sh" to run "genesis\_adaptive\_perm\_3chr.R"

## calculate empirical p-values

* "adaptive\_permutation\_summary.R" combines first 100 and adaptive
  permutations to calculate empirical p-values for original data (perm
  0) and 100 permutations

## calculate lambda GC (genomic inflation factor), principal components, and heritability

* calculate\_lambdaGC.R reads in all P values and calculates genomic inflation factor using a chisquare statistics
* GCTA\_heritability.sh makes a gRM in GCTA and uses that to calculate heritability in the original data and 100 permuted phenotype files
* PCA\_karyotypes.R calls inversion state and calculates principal components from hybrid swarm individuals

## identify index SNPs

* get\_clumping\_paramters.R reads through each set of permuted P values and identifies a P-value threshold to use in the plink clump command based on pulling out a top number of snps. The number of snps pulled out corresponds to the number of SNPs exceeding different FDR values in the original dataset
* plink\_clump\_perms.sh performs the index SNP calculation

## miscellaneous

*timeless\_genotypes.R extracts genotype information about the previously published indel in timeless using parental VCF (this polymorphism was not included in our GWAS)
