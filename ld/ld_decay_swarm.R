#!/usr/bin/env Rscript

#LD analysis
library(data.table)
library(SNPRelate)
library(foreach)
library(GWASTools)
#library(cowplot)
#library(viridis)
library(doMC)
registerDoMC(12)

#use imputed genos so missing data doesn't cause NA
print("hybrid swarm")
geno <- snpgdsOpen("/scratch/pae3g/genome-reconstruction/final2_draw1_replaced.vcf.gds", allow.fork=T)

a<-snpgdsSNPList(geno)
info<-data.table(snp.id=a$snp.id,
                 chr=a$chromosome,
                 pos=a$pos,
                 freq=a$afreq)
info[,maf:=pmin(freq, 1-freq)]
info<-info[maf>.1]

gwas<-fread("/scratch/pae3g/oldscratch_recovered/evolution/gwas_p_score_inv_id.txt")

gwas<-merge(gwas, info, by=c("chr", "pos"))

gwas[,fdr:=p.adjust(gwas.p, method="fdr")]

sample.id <- read.gdsn(index.gdsn(geno, "sample.id"))
swarm <- read.gdsn(index.gdsn(geno, "sample.annot/swarm"))
a.samps<-swarm=="A"
b.samps<-swarm=="B"
#look at long distance LD decay in random SNPs from swarm
y<-foreach(dist=c(50,100,500,1000,5000,10000,50000,100000,500000,1000000, 5000000,10000000), .errorhandling = "remove")%dopar%{   
    nsampled=1
    lddata<-list()
    while (nsampled<=10001){
        if(nsampled%%100==0){
            print(nsampled)
        }
        setkey(gwas, snp.id)
        focal.snp<-sample(gwas$snp.id, 1)
        #print(focal.snp)
        focal.chr<-gwas[J(focal.snp), chr]
        focal.pos<-gwas[J(focal.snp), pos]
        
        range<-gwas[chr==focal.chr&pos>(focal.pos+.95*dist)&pos<(focal.pos+1.05*dist)]
        #print(range)
        if (nrow(range)==0){
            next
        }
        candidates=range$snp.id
        if (nrow(range)==1) {
            focal.snp2=candidates
        } else{
            focal.snp2<-sample(candidates, 1)
        }
#this sampling was incorrect in previous iteration when there was only one possible snp
ld.a<-snpgdsLDMat(geno, slide=-1,verbose=F,snp.id = c(focal.snp, focal.snp2), sample.id=sample.id[a.samps], method='composite')$LD[2,1]
ld.b<-snpgdsLDMat(geno, slide=-1,verbose=F,snp.id = c(focal.snp, focal.snp2), sample.id=sample.id[b.samps], method='composite')$LD[2,1]
ld<-snpgdsLDMat(geno, slide=-1,verbose=F,snp.id = c(focal.snp, focal.snp2), method='composite')$LD[2,1]
lddata[[nsampled]]=data.table(draw=nsampled,
                              dist=dist,
                              ld.a=ld.a,
                              ld.b=ld.b,
                              ld=ld,
                              focal.snp=focal.snp,
                              focal.chr=focal.chr,
                              focal.snp2=focal.snp2,
                              focal.pos=focal.pos,
                              focal.pos2=range[snp.id==focal.snp2, pos])
nsampled = nsampled + 1

    }
    return(rbindlist(lddata))
}
y<-rbindlist(y)

write.table(y, "/scratch/pae3g/evolution/ld_decay_mafover0.1_byswarm_revised.txt", quote=F, sep="\t", row.names = F)

snpgdsClose(geno)

###############################################

##rare #####################################3

print("starting hybrid swarm rare")
geno <- snpgdsOpen("/scratch/pae3g/genome-reconstruction/final2_draw1_replaced.vcf.gds", allow.fork=T)

a<-snpgdsSNPList(geno)
info<-data.table(snp.id=a$snp.id,
                 chr=a$chromosome,
                 pos=a$pos,
                 freq=a$afreq)
info[,maf:=pmin(freq, 1-freq)]
info<-info[maf<.05]

gwas<-fread("/scratch/pae3g/oldscratch_recovered/evolution/gwas_p_score_inv_id.txt")

gwas<-merge(gwas, info, by=c("chr", "pos"))

gwas[,fdr:=p.adjust(gwas.p, method="fdr")]

sample.id <- read.gdsn(index.gdsn(geno, "sample.id"))
swarm <- read.gdsn(index.gdsn(geno, "sample.annot/swarm"))
a.samps<-swarm=="A"
b.samps<-swarm=="B"
#look at long distance LD decay in random SNPs from swarm
y<-foreach(dist=c(50,100,500,1000,5000,10000,50000,100000,500000,1000000, 5000000,10000000), .errorhandling = "remove")%dopar%{   
    nsampled=1
    lddata<-list()
    while (nsampled<=10001){
        if(nsampled%%100==0){
            print(nsampled)
        }
        setkey(gwas, snp.id)
        focal.snp<-sample(gwas$snp.id, 1)
        #print(focal.snp)
        focal.chr<-gwas[J(focal.snp), chr]
        focal.pos<-gwas[J(focal.snp), pos]
        
        range<-gwas[chr==focal.chr&pos>(focal.pos+.95*dist)&pos<(focal.pos+1.05*dist)]
        #print(range)
        if (nrow(range)==0){
            next
        }
        candidates=range$snp.id
        if (nrow(range)==1) {
            focal.snp2=candidates
        } else{
            focal.snp2<-sample(candidates, 1)
        }
        #this sampling was incorrect in previous iteration when there was only one possible snp
        ld.a<-snpgdsLDMat(geno, slide=-1,verbose=F,snp.id = c(focal.snp, focal.snp2), sample.id=sample.id[a.samps], method='composite')$LD[2,1]
        ld.b<-snpgdsLDMat(geno, slide=-1,verbose=F,snp.id = c(focal.snp, focal.snp2), sample.id=sample.id[b.samps], method='composite')$LD[2,1]
        ld<-snpgdsLDMat(geno, slide=-1,verbose=F,snp.id = c(focal.snp, focal.snp2), method='composite')$LD[2,1]
        lddata[[nsampled]]=data.table(draw=nsampled,
                                      dist=dist,
                                      ld.a=ld.a,
                                      ld.b=ld.b,
                                      ld=ld,
                                      focal.snp=focal.snp,
                                      focal.chr=focal.chr,
                                      focal.snp2=focal.snp2,
                                      focal.pos=focal.pos,
                                      focal.pos2=range[snp.id==focal.snp2, pos])
        nsampled = nsampled + 1
        
    }
    return(rbindlist(lddata))
}
y<-rbindlist(y)

write.table(y, "/scratch/pae3g/evolution/ld_decay_mafless0.05_byswarm_revised.txt", quote=F, sep="\t", row.names = F)

snpgdsClose(geno)

