#!/usr/bin/env Rscript

library(data.table)
library(SNPRelate)
library(gdsfmt)
library(doMC)
registerDoMC(20)
library(foreach)

phenos <- fread("/scratch/pae3g/oldscratch_recovered/phenos/phenos_062018.txt")
filters=fread("/scratch/pae3g/oldscratch_recovered/final_reconstruction2/hwe_missing_maf_filters.txt")
pass=filters[map_filter=="PASS", snp.id]

foreach(pop=c("A", "B"))%do%{
foreach(i=c(11:100))%do%{
    print(i)
    filestem<- paste0("/scratch/pae3g/genome-reconstruction/final2_draw", i, "_replaced")
    
    gds.fn=paste(filestem, ".vcf.gds", sep="")
    
    geno<-snpgdsOpen(gds.fn)
    #read genotypes and filters
    
    
    setkey(phenos, sample.id)
    
    genotyped_samples<- read.gdsn(index.gdsn(geno, "sample.id"))
    
    ids.to.use=intersect(genotyped_samples, phenos[swarm==pop]$sample.id)
    
    
    print("caclulating LOCO GRM")
    
    
    loco<-  snpgdsGRM(geno, method="Eigenstrat",sample.id=ids.to.use, snp.id=pass, autosome.only=FALSE, num.thread=20)
    rownames(loco$grm)=loco$sample.id
    colnames(loco$grm)=loco$sample.id
    
    a<-loco$grm
    save(a, file=paste(filestem,"_nonLOCO_Eigenstrat_allgwasSNPs_", pop, ".Rdat", sep=""))
    
    
    snpgdsClose(geno)
}
}

