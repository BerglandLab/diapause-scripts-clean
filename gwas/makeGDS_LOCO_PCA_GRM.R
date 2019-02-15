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


    pcaOut <- foreach(i=c("chr2L", "chr2R", "chr3L" ,"chr3R", "chrX"))%do%{
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
    pca.temp <- snpgdsPCA(geno, snp.id=snpset.use, autosome.only=FALSE, maf=maf, sample.id = ids.to.use,
                          missing.rate=missing.rate, num.thread=threads)
    
  #snpgdsPCASNPLoading(pca.temp, genofile)
  
    pca.dt <- data.table(pc1=pca.temp$eigenvect[,1],
                       pc2=pca.temp$eigenvect[,2],
                       pc3=pca.temp$eigenvect[,3],
                       pc4=pca.temp$eigenvect[,4],
                       pc5=pca.temp$eigenvect[,5],
                       pc6=pca.temp$eigenvect[,6],
                       pc7=pca.temp$eigenvect[,7],
                       pc8=pca.temp$eigenvect[,8],
                       pc9=pca.temp$eigenvect[,9],
                       pc10=pca.temp$eigenvect[,10],
                       pc11=pca.temp$eigenvect[,11],
                       pc12=pca.temp$eigenvect[,12],
                       pc13=pca.temp$eigenvect[,13],
                       pc14=pca.temp$eigenvect[,14],
                       pc15=pca.temp$eigenvect[,15],
                       pc16=pca.temp$eigenvect[,16],
                       pc17=pca.temp$eigenvect[,17],
                       pc18=pca.temp$eigenvect[,18],
                       pc19=pca.temp$eigenvect[,19],
                       pc20=pca.temp$eigenvect[,20],
                       pc21=pca.temp$eigenvect[,21],
                       pc22=pca.temp$eigenvect[,22],
                       pc23=pca.temp$eigenvect[,23],
                       pc24=pca.temp$eigenvect[,24],
                       pc25=pca.temp$eigenvect[,25],
                       pc26=pca.temp$eigenvect[,26],
                       pc27=pca.temp$eigenvect[,27],
                       pc28=pca.temp$eigenvect[,28],
                       pc29=pca.temp$eigenvect[,29],
                       pc30=pca.temp$eigenvect[,30],
                       pc31=pca.temp$eigenvect[,31],
                       pc32=pca.temp$eigenvect[,32],
                       id=pca.temp$sample,
                       chrs=i)
    #pca.dt
  }
  
  pcaOut <- rbindlist(pcaOut)

  write.table(pcaOut, paste(filestem, "_LOCO_PCA_strict_both.txt", sep=""), quote=FALSE, row.names=FALSE, sep="\t")

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


