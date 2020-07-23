
#!/usr/bin/env Rscript
library(data.table)
library(SNPRelate)
library(foreach)
library(GWASTools)
library(doMC)
registerDoMC(20)

#use imputed genos so missing data doesn't cause NA
phenos <- fread("/scratch/pae3g/revisions/phenos_wolbachia_missinggeno.txt")
#filter missing data
phenos <- phenos[n.chr.imp==5&prop.unknown<=0.05]


foreach(pop=c("A", "B", "both")) %do%{
    #define samples
    if(pop=="both"){
        samps=phenos$sample.id
    }else{
        samps=phenos[swarm==pop, sample.id]
    }
    geno <- snpgdsOpen("/scratch/pae3g/genome-reconstruction/final2_draw1_replaced.vcf.gds", allow.fork=T)
    
    #get snp info for focal population  
    a<-snpgdsSNPList(geno, sample.id=samps)
    info<-data.table(snp.id=a$snp.id,
                     chr=a$chromosome,
                     pos=a$pos,
                     freq=a$afreq)
    info[,maf:=pmin(freq, 1-freq)]
    filters=fread("/scratch/pae3g/oldscratch_recovered/final_reconstruction2/hwe_missing_maf_filters.txt")
    
    info<-merge(info, filters[,.(chr, pos, map_filter)], by=c("chr", "pos"))
    snpstouse<-info[maf>=0.05&map_filter=="PASS", .(chr, pos)]
    setkey(snpstouse, chr, pos)
    snpgdsClose(geno)
    
    
    geno <- snpgdsOpen("/scratch/pae3g/revisions/hs.hc.snp.99.9.both.linenames.vcf.gds", allow.fork=T)
    
    
    a<-snpgdsSNPList(geno)
    info<-data.table(snp.id=a$snp.id,
                     chr=a$chromosome,
                     pos=a$pos)
    setkey(info, chr, pos)
    info<-merge(info, snpstouse)
    sample.id <- read.gdsn(index.gdsn(geno, "sample.id"))
    parents<-fread("/scratch/pae3g/overwintering/parent_line_by_swarm.txt", header=F)
    setnames(parents, c("line","pid1", "pid2", "swarm"))
    
    #filter out single swarm samples
    if(pop!="both"){
        samplestouse<-parents[swarm==pop, line]
    }else{
        samplestouse<-parents[,line]
    }
    
    m<-foreach(draw=c(1:2), .errorhandling="remove")%dopar%{
        if(draw%%100==0){
            print(draw)
        }
        #print(draw)
        y<-foreach(chrom=c(2,3))%do%{
            setkey(info, snp.id)
            focal.snp<-sample(info[chr==paste0(chrom,"L"), snp.id], 1)
            #print(focal.snp)
            match.snp<-sample(info[chr==paste0(chrom,"R"), snp.id], 1)
            ld<-snpgdsLDMat(geno, slide=-1,verbose=F,snp.id = c(focal.snp, match.snp),sample.id=samplestouse, method='composite')$LD[2,1]
            return(data.table(draw=draw,
                              chr=chrom,
                              ld=ld,
                              focal.snp=focal.snp,
                              match.snp=match.snp))
        }
        return(rbindlist(y))
    }
    m<-rbindlist(m)
    
    write.table(m, paste0("/scratch/pae3g/revisions/chromsome_arm_ld_parents_pop", pop, ".txt"), quote=F, sep="\t", row.names=F)
    
    
    m<-foreach(draw=c(1:2), .errorhandling="remove")%dopar%{
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
            ld<-snpgdsLDMat(geno, slide=-1,verbose=F,snp.id = c(focal.snp, match.snp), sample.id=samplestouse, method='composite')$LD[2,1]
            return(data.table(draw=draw,
                              chr=chrom,
                              ld=ld,
                              focal.snp=focal.snp,
                              match.snp=match.snp))
        }
        return(rbindlist(y))
    }
    m<-rbindlist(m)
    
    write.table(m, paste0("/scratch/pae3g/revisions/interchromosomal_ld_parents_pop", pop, ".txt"), quote=F, sep="\t", row.names=F)
    snpgdsClose(geno)
}

