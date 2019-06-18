# Overview

These scripts perform the population genetic analysis shown in Figure 8

## IHS
* ihs.sh preps the DGRP file to use in the rehh package in R
* ihs\_dgrp.R calculates ihs genome wide and then pulls out ihs for index SNPs in all permutations

## Tajima's D
This code is modified from Ruzicka et al 2019 (PLOS Biology)

* tajimas\_d.sh prepares the dgrp vcf file (available from teh dgrp website) to use in the PopGenome package in R.
* tajimas\_d\_dgrp.R calculates sliding window Tajima's D and then pulls out windows with index SNPs for the analysis

## DPGP3
This script uses a file with all DPGP allele frequencies that were extracted from teh fasta files via a custom script

* dpgp\_zambia.R extracts pro-diapause allele frequencies for all index SNPs
