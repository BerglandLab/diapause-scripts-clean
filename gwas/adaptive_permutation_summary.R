
library(data.table)
library(foreach)
library(cowplot)
library(doMC)
registerDoMC(16)

#read in snp metadata
dat=fread("/scratch/pae3g/final_reconstruction2/hwe_missing_maf_filters.txt")

#read in summarized permutation data from teh initial 100 permutations of 100 imputations
initial_perms<-foreach(i=c(101:200)) %dopar% {
  x<-fread(paste("/scratch/pae3g/final_reconstruction2/perm", i,"_100_draw_summary.txt", sep=""))
  x[,perm:=i]
  return(x)
}

initial_perms<-rbindlist(initial_perms)

#read in the summarized permutation data from 500 adaptive permutations of the less common snp classes
adaptive_perms<-foreach(i=c(1001:1500), .errorhandling = "remove") %dopar% {
  x<-fread(paste("/scratch/pae3g/final_reconstruction2/adaptive_perm_", i,".txt", sep=""))
  x[,perm:=i]
  return(x)
}

adaptive_perms<-rbindlist(adaptive_perms)

#get rid of max.pval for the merge because it is only in one dataset
adaptive_perms[, max.pval:=NULL]

#make names match
setnames(adaptive_perms, "avg.pval", "mean.p")
setnames(adaptive_perms, "min.pval", "min.p")
setnames(adaptive_perms, "med.pval", "med.p")
setnames(adaptive_perms, "sd.pval", "sd.p")

setnames(initial_perms, "sd.pval", "sd.p")


#merge in seeds of adaptive_perms
input=fread("/scratch/pae3g/genome-reconstruction/seeded_permutation_input.txt")
names(input)=c("file", "pop", "pheno", "draw", "perm", "seed")

initial_perms<-merge(initial_perms, input[draw==1,.(perm, seed)], by="perm")

#rbind the two data sets

perms<-rbind(adaptive_perms, initial_perms)

rm(adaptive_perms)
rm(initial_perms)

#merge permutation data with metadata

setkey(perms, snpID)
setkey(dat, snpID)

perms<-merge(perms, dat)

#read in actual data

actual.sum=fread("/scratch/pae3g/final_reconstruction2/100draws_actual_sum.csv")

actual.sum=merge(actual.sum, dat, by=c("snpID", "chr", "pos"))


#calculate percentile-ranked pvalues compared to same-bin snps
pvals<-foreach(maf=unique(actual.sum$maf.bin))%do%{
  print(maf)
  missing.bins<-foreach(mis=unique(actual.sum$mis.bin), .errorhandling="remove")%do%{
    print(mis)
    current.snps=actual.sum[maf.bin==maf&mis.bin==mis]
    current.ecdf=ecdf(perms[maf.bin==maf&mis.bin==mis, mean.p])
    current.snps[,pval.percentile:=current.ecdf(avg.pval)]
    current.snps[,n.permuted.snps:=dim(perms[maf.bin==maf&mis.bin==mis])[1]]
    return(current.snps)
  }
  missing.bins<-rbindlist(missing.bins)
  return(missing.bins)
}

pvals<-rbindlist(pvals)
pvals[,V1:=NULL]
pvals[,lod:=-log10(avg.pval)]
#replace zero pvals with 1/n 
pvals[pval.percentile==0, pval.percentile:=1/n.permuted.snps]
pvals[,lod.percentile:=-log10(pval.percentile)]
pvals[,fdr.percentile:=p.adjust(pval.percentile, method="fdr")]

#final data set
write.table(pvals[,.(chr, pos, snpID, avg.pval, maf.bin, mis.bin, MAF, lod, pval.percentile, lod.percentile, fdr.percentile)], "/scratch/pae3g/final_reconstruction2/adaptive_perms_pvalues.txt", sep="\t", quote=F)

#now do same calculation of a percentile-ranked pvalue for each permutation and calculate lambdaGC for each permutation

lambda.table<-foreach(perm.num=c(101:200), .combine="rbind")%do%{ 
  print(perm.num)
  #pull out that permutation data
  actual.sum=perms[perm==perm.num]
  pvals<-foreach(maf=unique(actual.sum$maf.bin))%do%{
    print(maf)
    missing.bins<-foreach(mis=unique(actual.sum$mis.bin), .errorhandling="remove")%do%{
      print(mis)
      #remove the current permutation
      current.snps=actual.sum[maf.bin==maf&mis.bin==mis]
      current.ecdf=ecdf(perms[perm!=perm.num&maf.bin==maf&mis.bin==mis, mean.p])
      current.snps[,pval.percentile:=current.ecdf(mean.p)]
      current.snps[,n.permuted.snps:=dim(perms[maf.bin==maf&mis.bin==mis])[1]]
      return(current.snps)
    }
    missing.bins<-rbindlist(missing.bins)
    return(missing.bins)
  }

  pvals<-rbindlist(pvals)
  pvals[,lod:=-log10(mean.p)]
  pvals[,lod.percentile:=-log10(pval.percentile)]
  pvals[,percentile.adjust:=p.adjust(pval.percentile, method="fdr")]

  write.table(pvals, paste("/scratch/pae3g/final_reconstruction2/adaptive_perms_perm", perm.num, "_pvalues.txt",sep=""), sep="\t", quote=F)
  
  chisq<-qchisq(1-pvals$pval.percentile,1)
  
  return(data.table(lambda=median(chisq)/qchisq(0.5,1), perm=perm.num))
}

#did not get full table because connection was interrupted, will have to read in and caculate again
write.csv(lambda.table, "/scratch/pae3g/final_reconstruction2/perm101-200lambdaGC.csv")

perm.percentile<-foreach(perm.num=c(101:200))%dopar%{
  return(fread(paste("/scratch/pae3g/final_reconstruction2/adaptive_perms_perm", perm.num, "_pvalues.txt",sep="")))
}
perm.percentile<-rbindlist(perm.percentile)
lambda<-perm.percentile[,.(lambda=median(qchisq(1-pval.percentile,1))/qchisq(0.5,1)), .(perm)]

actual.lambda<-median(qchisq(1-pvals$pval.percentile,1))/qchisq(0.5,1)

ggplot(perm.percentile, aes(x=pval.percentile))+geom_histogram(bins = 1000)+facet_wrap(~perm)
library(data.table)
library(cowplot)

#analyze the results of the adaptive permutation (100 initial permutations plus 500 permutations of snps in rarer categories)
pvals<-fread("/scratch/pae3g/final_reconstruction2/adaptive_perms_pvalues.txt")

#histogram of pvalues
ggplot(pvals, aes(x=pval.percentile))+geom_histogram(bins = 1000)




