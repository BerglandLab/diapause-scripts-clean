#LD analysis
library(data.table)
library(SNPRelate)
library(foreach)
library(GWASTools)
#library(cowplot)
#library(viridis)
library(doMC)
registerDoMC(12)


#############################
#compare to dgrp gds'
#############################


print("starting dgrp common")

#
geno<-snpgdsOpen("/scratch/pae3g/oldscratch_recovered/evolution/dgrp2.vcf.gds", allow.fork=T)
a<-snpgdsSNPList(gdsobj=geno)
samps<-snpgdsSummary(geno)$sample.id

info<-data.table(snp.id=a$snp.id,
                 chr=a$chromosome,
                 pos=a$pos,
                 freq=a$afreq)
info[,maf:=pmin(freq, 1-freq)]
info<-info[maf>.1]

y<-foreach(dist=c(50,100,500,1000,5000,10000,50000,100000,500000,1000000,5000000,10000000), .errorhandling = "remove")%dopar%{   
    nsampled=1
    lddata<-list()
    while (nsampled<=10000){
        
        #while (nsampled<=10000){
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
        
        ld<-snpgdsLDMat(geno, slide=-1,verbose=F,snp.id = c(focal.snp, focal.snp2), sample.id=sample(samps, 68), method='composite')$LD[2,1]
        lddata[[nsampled]]=data.table(draw=nsampled,
                                      dist=dist,
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

write.table(y, "/scratch/pae3g/evolution/ld_decay_dgrp_mafover0.1_downsampled_revised.txt", quote=F, sep="\t", row.names = F)

snpgdsClose(geno)


########################################
#dgrp rare

print("starting dgrp rare")
geno<-snpgdsOpen("/scratch/pae3g/oldscratch_recovered/evolution/dgrp2.vcf.gds", allow.fork=T)
a<-snpgdsSNPList(gdsobj=geno)
samps<-snpgdsSummary(geno)$sample.id

info<-data.table(snp.id=a$snp.id,
                 chr=a$chromosome,
                 pos=a$pos,
                 freq=a$afreq)
info[,maf:=pmin(freq, 1-freq)]
info<-info[maf<.05]

y<-foreach(dist=c(50,100,500,1000,5000,10000,50000,100000,500000,1000000, 5000000,10000000), .errorhandling = "remove")%dopar%{   
    nsampled=1
    lddata<-list()
    while (nsampled<=10000){
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
        
        ld<-snpgdsLDMat(geno, slide=-1,verbose=F,snp.id = c(focal.snp, focal.snp2), sample.id=sample(samps, 68), method='composite')$LD[2,1]
        lddata[[nsampled]]=data.table(draw=nsampled,
                                      dist=dist,
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

write.table(y, "/scratch/pae3g/evolution/ld_decay_dgrp_mafless0.05_downsampled_revised.txt", quote=F, sep="\t", row.names = F)