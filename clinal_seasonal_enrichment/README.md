# Overview

This folder contains the scripts used to perform the clinal and seasonal enrichment tests and polygenic score calculations using data from Bergland et al 2014 and Machado et al 2018. The clinal and seasonal P values are available on Data Dryad. All tests were performed using index SNPs, also available on Dryad.

## enrichment tests

* enrich\_bergland2014.R
* enrich\_machado2018.R
* enrich\_machado2018\_bypop.R

## test for parallel SNP usage in individual populations

* three\_population\_overlap\_chisq.R

## polygenic score tests (multiplying effect sizes)

* polygenic\_scores.R
* polygenic\_score\_by\_population.R

## Stouffer's combined P value test for GWAS and seasonal data
* stouffers\_combinedP.R
