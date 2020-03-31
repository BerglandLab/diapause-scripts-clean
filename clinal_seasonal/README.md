# Overview

This folder contains the scripts used to perform  polygenic score calculations using data from Bergland et al 2014 and Machado et al 2019.  All tests were performed using LASSO SNPs as well as threshold-ranked SNPs

## polygenic score tests (multiplying effect sizes)

*core20delta.R calculates the individual population changes in allele frequencyes (spring-fal, logit transformed)

* sign\_test\_universal\_threshold.R intersects the quantile-ranked GWAS data with two clinal adn seasonal datasets to calculate polygenic scores

* analysis on LASSO SNPs is done in the figures-generating R script
* clinal and seasonal P values available on Dryad
