library(rehh)
library(data.table)
library(foreach)
library(doMC)
registerDoMC(20)
#fix input files
foreach(chr=c("2L", "2R", "3L", "3R", "X"))%do%{
     inp<-fread(paste0("/mnt/pricey_2/priscilla/dgrp2.filtered.", chr, ".impute.legend"), header=T)
     inp[,chr:=tstrsplit(ID, split="_")[[1]]]
     inp[,allele0:=0]
     inp[,allele1:=1]
     write.table(inp[,.(ID, chr, pos, allele0, allele1)], paste0("/mnt/pricey_2/priscilla/dgrp2.", chr, ".inp"), quote=F, sep=" ", row.names=F, col.names=F)
 }
 
 a<-foreach(chr=c("2L", "2R", "3L", "3R", "X"))%do%{
     hap<-data2haplohh(hap_file=paste0("dgrp2.filtered.", chr,".impute.hap"),map_file=paste0("dgrp2.", chr, ".inp"),min_perc_geno.hap=90, min_perc_geno.snp=90, recode.allele=TRUE, haplotype.in.columns=TRUE)
 
     res<-scan_hh(hap, threads=20)
     return(as.data.table(res))
 }
 a<-rbindlist(a)
 
 write.table(a, "/mnt/pricey_2/priscilla/ihs.txt", quote=F, sep="\t", row.names=F)

#move ihs.txt to rivanna to work there

