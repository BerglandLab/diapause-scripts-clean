# Overview

These scripts perform the population genetic analysis shown in Figures 7 and 8

## Northern lines

* The folder "northern" contains 5 scripts that were used to map the reads downloaded from the SRA then call SNPs in gVCFs, build a VCF, and filter that VCF

## IHS
ihs was analyzed relative to two north american populations: the DGRP (north clrolina) and a set of 205 lines collected in Pennsylvania and Maine ("northern"). All data on dryad
* "ihs\_dgrp.sh" and "ihs\_northern.sh" preps the file sto use in the rehh package in R
* "ihs\_dgrp.R" and "ihs\_northern.R" calculates the haplotype files
* "ihs\_universal\_lasso.R" calculates ihs for all LASSO SNPs
* "ihs\_universal\_threshold.R" calculates ihs for quantile ranked SNPs from the gwas

## DPGP3

This analysis uses a file with all DPGP allele frequencies that were extracted from teh fasta files via a custom script https://github.com/alanbergland/DEST/

## African admixture

This analysis tested whether diapause-associated SNPs were found in regions of European admixture in Zambian flies
* "admix\_universal\_lasso.R" caculates admixture for LASSO SNPs
* "admix\_univerals\_threshold.R" calculates admixture for quantile-ranked thresholds
