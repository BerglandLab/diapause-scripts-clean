
#!/usr/bin/env Rscript
library(data.table)
library(SNPRelate)
library(foreach)
library(GWASTools)
library(doMC)
registerDoMC(20)

#use imputed genos so missing data doesn't cause NA
print("dgrp")
geno <- snpgdsOpen("/scratch/pae3g/oldscratch_recovered/evolution/dgrp2.vcf.gds", allow.fork=T)

a<-snpgdsSNPList(geno)
info<-data.table(snp.id=a$snp.id,
                 chr=a$chromosome,
                 pos=a$pos,
                 freq=a$afreq)
info[,maf:=pmin(freq, 1-freq)]

#test 10,000 pairs for within-chromosome long-range LD (2L to 2R; 3L to 3R)

m<-foreach(draw=c(1:10000), .errorhandling="remove")%dopar%{
    if(draw%%100==0){
        print(draw)
    }
    #print(draw)
    y<-foreach(chrom=c(2,3))%do%{
        setkey(info, snp.id)
        focal.snp<-sample(info[chr==paste0(chrom,"L"), snp.id], 1)
        #print(focal.snp)
        match.snp<-sample(info[chr==paste0(chrom,"R"), snp.id], 1)
        ld<-snpgdsLDMat(geno, slide=-1,verbose=F,snp.id = c(focal.snp, match.snp), method='composite')$LD[2,1]
        return(data.table(draw=draw,
                          chr=chrom,
                          ld=ld,
                          focal.snp=focal.snp,
                          match.snp=match.snp))
    }
    return(rbindlist(y))
}
m<-rbindlist(m)

write.table(m, "/scratch/pae3g/revisions/chromsome_arm_ld_dgrp.txt", quote=F, sep="\t", row.names=F)


m<-foreach(draw=c(1:10000), .errorhandling="remove")%dopar%{
    if(draw%%100==0){
        print(draw)
    }
    #print(draw)
    y<-foreach(chrom=c(2,3, "X"))%do%{
        setkey(info, snp.id)
        if(chrom=="2"){
            focal.snp<-sample(info[chr==2L|chr=="2R", snp.id], 1)
        #print(focal.snp)
            match.snp<-sample(info[chr=="3L"|chr=="3R"|chr=="X", snp.id], 1)
            
        }
        else if(chrom=="3"){
            focal.snp<-sample(info[chr==3L|chr=="3R", snp.id], 1)
            #print(focal.snp)
            match.snp<-sample(info[chr=="2L"|chr=="2R"|chr=="X", snp.id], 1)
            
        }
        else if(chrom=="X"){
            focal.snp<-sample(info[chr=="X", snp.id], 1)
            #print(focal.snp)
            match.snp<-sample(info[chr=="2L"|chr=="2R"|chr=="3L"|chr=="3R", snp.id], 1)
        }
        ld<-snpgdsLDMat(geno, slide=-1,verbose=F,snp.id = c(focal.snp, match.snp), method='composite')$LD[2,1]
        return(data.table(draw=draw,
                          chr=chrom,
                          ld=ld,
                          focal.snp=focal.snp,
                          match.snp=match.snp))
    }
    return(rbindlist(y))
}
m<-rbindlist(m)

write.table(m, "/scratch/pae3g/revisions/interchromosomal_ld_dgrp.txt", quote=F, sep="\t", row.names=F)

