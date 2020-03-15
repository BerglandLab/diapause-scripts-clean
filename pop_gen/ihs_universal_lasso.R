# library(rehh)
# library(data.table)
# library(foreach)
# library(doMC)
# registerDoMC(20)
#fix input files
# foreach(chr=c("2L", "2R", "3L", "3R", "X"))%do%{
#      inp<-fread(paste0("/mnt/pricey_2/priscilla/dgrp2.filtered.", chr, ".impute.legend"), header=T)
#      inp[,chr:=tstrsplit(ID, split="_")[[1]]]
#      inp[,allele0:=0]
#      inp[,allele1:=1]
#      write.table(inp[,.(ID, chr, pos, allele0, allele1)], paste0("/mnt/pricey_2/priscilla/dgrp2.", chr, ".inp"), quote=F, sep=" ", row.names=F, col.names=F)
#  }
#  
#  a<-foreach(chr=c("2L", "2R", "3L", "3R", "X"))%do%{
#      hap<-data2haplohh(hap_file=paste0("dgrp2.filtered.", chr,".impute.hap"),map_file=paste0("dgrp2.", chr, ".inp"),min_perc_geno.hap=90, min_perc_geno.snp=90, recode.allele=TRUE, haplotype.in.columns=TRUE)
#  
#      res<-scan_hh(hap, threads=20)
#      return(as.data.table(res))
#  }
#  a<-rbindlist(a)
#  
#  write.table(a, "/mnt/pricey_2/priscilla/ihs.txt", quote=F, sep="\t", row.names=F)

#move ihs.txt to rivanna to work there

library(data.table)
library(foreach)
library(rehh)
library(doMC)
registerDoMC(20)

#read ihs. data
a.raw<-fread("/scratch/pae3g/revisions/evolution/ihs.txt")

#read file info to loop through
files<-fread("/scratch/pae3g/genome-reconstruction/universal_input.txt")

#read in files
p<-foreach(pop=files$V1[1:3000], phenotype=files$V2[1:3000],draw=files$V3[1:3000], perm=files$V4[1:3000], model=files$V6[1:3000], .errorhandling="remove")%dopar%{
    load(paste0("/scratch/pae3g/revisions/lasso/", "lasso_", phenotype, "_draw",draw, "_perm", perm, "_", model, "_pop", pop, ".lassoSites.Rdata"))
    sites[, chr:=tstrsplit(site,split="_")[[1]]]
    sites[, pos:=as.numeric(tstrsplit(site,split="_")[[2]])]
    print(paste(pop, phenotype, draw, perm, model, sep=", "))
    
    #read in GWAS p values and scores
    load(paste("/scratch/pae3g/revisions/genesis_", phenotype, "_draw", draw, "_perm", perm, "_pop", pop, "_" , model, "_allsnpgrm_wolbachia.Rdat", sep=""))

    #merge wtih info because chromosomes got lost in genesis
    gwas<-assoc.results
    setnames(gwas, c("chr", "pos"), c("CHR", "POSITION"))

    #merge with ihs data
    a<-merge(a.raw, gwas, by=c("CHR", "POSITION"), all.x=T)
    #currently coded so that ancestral allele is reference (0) and derived allele is alternate (1)
    #we want "derived" to be the pro-diapause allele
    #positive score means that ref is pro-diapause (aka ancestral is pro diapause)
    #so, when gwas score is POSTIVE, switch ihh_a and ihh_d before computing ihs

    a[Score>=0,ihh.a:=iHH_D]
    a[is.na(Score)|Score<0,ihh.a:=iHH_A]
    a[Score>=0,ihh.d:=iHH_A]
    a[is.na(Score)|Score<0,ihh.d:=iHH_D]

    #flip allele frequency if ancestral and derived are flipped
    a[Score>=0, freq_A:=1-freq_A]

    #delete original calls and rename
    a[,iHH_A:=NULL]
    a[,iHH_D:=NULL]
    setnames(a, c("ihh.a", "ihh.d"), c("iHH_A", "iHH_D"))
    a<-a[,.(CHR, POSITION, freq_A,  iHH_A, iHH_D, iES_Tang_et_al_2007, iES_Sabeti_et_al_2007)]

    b<-ihh2ihs(a) #note the default is discard maf < 0.05
    c<-as.data.table(b$ihs)

    d<-merge(gwas, c, by=c("CHR", "POSITION"))
    setnames(d, c("CHR","POSITION", "LOGPVALUE"), c("chr", "pos",'ihs.p'))

    f<-merge(sites[,.(chr, pos)], d, by=c("chr", "pos"))
    #pull out most significant snp in each block
    f[,perm:=perm]
    f[,draw:=draw]
    f[,pheno:=phenotype]
    f[,pop:=pop]
    f[,model:=model]
    f[,test:="dgrp"]
    #save output
    return(f)
}

p<-rbindlist(p)


save(p, file="/scratch/pae3g/revisions/evolution/ihs_dgrp_universal_lasso.Rdata")


a.raw<-fread("/scratch/pae3g/revisions/evolution/BME_LN_LNPA.ihs.txt")

p<-foreach(pop=files$V1[1:3000], phenotype=files$V2[1:3000],draw=files$V3[1:3000], perm=files$V4[1:3000], model=files$V6[1:3000], .errorhandling="remove")%dopar%{
    load(paste0("/scratch/pae3g/revisions/lasso/", "lasso_", phenotype, "_draw",draw, "_perm", perm, "_", model, "_pop", pop, ".lassoSites.Rdata"))
    sites[, chr:=tstrsplit(site,split="_")[[1]]]
    sites[, pos:=as.numeric(tstrsplit(site,split="_")[[2]])]
    print(paste(pop, phenotype, draw, perm, model, sep=", "))
    
    #read in GWAS p values and scores
    load(paste("/scratch/pae3g/revisions/genesis_", phenotype, "_draw", draw, "_perm", perm, "_pop", pop, "_" , model, "_allsnpgrm_wolbachia.Rdat", sep=""))
    
    #merge wtih info because chromosomes got lost in genesis
    gwas<-assoc.results
    setnames(gwas, c("chr", "pos"), c("CHR", "POSITION"))
    
    #merge with ihs data
    a<-merge(a.raw, gwas, by=c("CHR", "POSITION"), all.x=T)
    #currently coded so that ancestral allele is reference (0) and derived allele is alternate (1)
    #we want "derived" to be the pro-diapause allele
    #positive score means that ref is pro-diapause (aka ancestral is pro diapause)
    #so, when gwas score is POSTIVE, switch ihh_a and ihh_d before computing ihs
    
    a[Score>=0,ihh.a:=iHH_D]
    a[is.na(Score)|Score<0,ihh.a:=iHH_A]
    a[Score>=0,ihh.d:=iHH_A]
    a[is.na(Score)|Score<0,ihh.d:=iHH_D]
    
    #flip allele frequency if ancestral and derived are flipped
    a[Score>=0, freq_A:=1-freq_A]
    
    #delete original calls and rename
    a[,iHH_A:=NULL]
    a[,iHH_D:=NULL]
    setnames(a, c("ihh.a", "ihh.d"), c("iHH_A", "iHH_D"))
    a<-a[,.(CHR, POSITION, freq_A,  iHH_A, iHH_D, iES_Tang_et_al_2007, iES_Sabeti_et_al_2007)]
    
    b<-ihh2ihs(a) #note the default is discard maf < 0.05
    c<-as.data.table(b$ihs)
    
    d<-merge(gwas, c, by=c("CHR", "POSITION"))
    setnames(d, c("CHR","POSITION", "LOGPVALUE"), c("chr", "pos",'ihs.p'))
    
    f<-merge(sites[,.(chr, pos)], d, by=c("chr", "pos"))
    #pull out most significant snp in each block
    f[,perm:=perm]
    f[,draw:=draw]
    f[,pheno:=phenotype]
    f[,pop:=pop]
    f[,model:=model]
    f[,test:="north"]
    #save output
    return(f)
}

p<-rbindlist(p)


save(p, file="/scratch/pae3g/revisions/evolution/ihs_northern_universal_lasso.Rdata")
