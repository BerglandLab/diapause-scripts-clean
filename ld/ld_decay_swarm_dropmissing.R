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
phenos <- fread("/scratch/pae3g/revisions/phenos_wolbachia_missinggeno.txt")
#filter missing data
phenos <- phenos[n.chr.imp==5&prop.unknown<=0.05]


foreach(pop=c("A", "B", "both")) %do% {
    #define samples
if(pop=="both"){
    samps=phenos$sample.id
}else{
    samps=phenos[swarm==pop, sample.id]
}
 
#get snp info for focal population  
a<-snpgdsSNPList(geno, sample.id=samps)
info<-data.table(snp.id=a$snp.id,
                 chr=a$chromosome,
                 pos=a$pos,
                 freq=a$afreq)
info[,maf:=pmin(freq, 1-freq)]
filters=fread("/scratch/pae3g/oldscratch_recovered/final_reconstruction2/hwe_missing_maf_filters.txt")

info<-merge(info, filters[,.(chr, pos, map_filter)], by=c("chr", "pos"))
snps<-info[maf>=0.05&map_filter=="PASS", snp.id]
#rm(filters)

#filter out individuals wtih missing reconstrucion data as in GWAS

#look at long distance LD decay in random SNPs from swarm
y<-foreach(dist=c(50,100,500,1000,5000,10000,50000,100000,500000,1000000, 5000000,10000000), .errorhandling = "remove")%dopar%{   
    nsampled=1
    lddata<-list()
    while (nsampled<=10001){
        if(nsampled%%100==0){
            print(nsampled)
        }
        setkey(info, snp.id)
        focal.snp<-sample(info$snp.id, 1)
        #print(focal.snp)
        focal.chr<-info[J(focal.snp), chr]
        focal.pos<-info[J(focal.snp), pos]
        
        range<-info[chr==focal.chr&pos>(focal.pos+.95*dist)&pos<(focal.pos+1.05*dist)]
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

ld<-snpgdsLDMat(geno, slide=-1,verbose=F,sample.id=samps, snp.id = c(focal.snp, focal.snp2), method='composite')$LD[2,1]
lddata[[nsampled]]=data.table(draw=nsampled,
                              dist=dist,
                              pop=pop,
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

write.table(y, paste0("/scratch/pae3g/evolution/ld_decay_maf0.05_pop", pop, ".txt"), quote=F, sep="\t", row.names = F)

}
