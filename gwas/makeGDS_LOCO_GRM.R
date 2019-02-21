#!/usr/bin/env Rscript

library(GWASTools)
library(data.table)
library(SNPRelate)
library(gdsfmt)
library(doMC)
registerDoMC(20)
library(foreach)


args=commandArgs(trailingOnly=TRUE)

filestem<- args[1]

print("making GDS!")
vcf.fn=paste(filestem, ".vcf", sep="")
gds.fn=paste(vcf.fn, ".gds", sep="")

#check if a gds has already been made and make one if not.
if (!file.exists(gds.fn)) {
  snpgdsVCF2GDS(vcf.fn, gds.fn, method=c("biallelic.only"), snpfirstdim = FALSE)
} else{
    print("GDS already exists, continuing processing")
  }


print("opening GDS!")
geno<-snpgdsOpen(gds.fn, allow.fork=TRUE, readonly = FALSE)

#read genotypes and filters
print("filtering snps")
filters=fread("/scratch/pae3g/final_reconstruction2/hwe_missing_maf_filters.txt")


maf <- 0.05
missing.rate <- 0.15
threads <- 20
pass=filters[qc_filter=="PASS", snp.id]

phenos <- fread("/scratch/pae3g/phenos/phenos_062018.txt")
setkey(phenos, sample.id)

genotyped_samples<- read.gdsn(index.gdsn(geno, "sample.id"))
  #only use genotypes from swarm

if(!file.exists(paste(filestem,"_LOCO_Eigenstrat_strict_both.Rdat", sep="")) | !file.exists(paste(filestem, "_LOCO_PCA_strict_both.txt", sep=""))) {
    
    print("caculating PCA")
    ids.to.use=intersect(genotyped_samples, phenos$sample.id)

    snpset <- snpgdsLDpruning(geno, 
                          ld.threshold=0.2, 
                          slide.max.bp = 5000, 
                          autosome.only=FALSE,
                          missing.rate=.15,
                          sample.id=ids.to.use,
                          snp.id=pass,
                          maf=.05)

  print("caclulating LOCO GRM")

  loco.list=list()
  foreach(i=c( "chr2L", "chr2R", "chr3L" ,"chr3R", "chrX"))%do%{
    print(i)
    if(i=="chr2L" | i=="chr2R"){
      chr.to.use=c("chr3L", "chr3R", "chrX")
    }
    else if(i=="chr3L"| i=="chr3R"){
      chr.to.use=c("chr2L", "chr2R", "chrX")
    }
    else if(i=="chrX"){
      chr.to.use=c("chr2L", "chr2R", "chr3L", "chr3R")
    }
    snpset.use <- unlist(snpset[chr.to.use])
    print(length(snpset.use))
    loco<-  snpgdsGRM(geno, method="Eigenstrat", snp.id=snpset.use, sample.id=ids.to.use, autosome.only=FALSE, num.thread=20)
    loco.list[[i]]=loco$grm
    rownames(loco.list[[i]])=loco$sample.id
    colnames(loco.list[[i]])=loco$sample.id
  }
  
  save(loco.list, file=paste(filestem,"_LOCO_Eigenstrat_strict_both.Rdat", sep=""))
  } else{
    print("LOCO PCA and GRM already computed for this file")
  }

snpgdsClose(geno)


