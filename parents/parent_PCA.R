library(data.table)
library(SNPRelate)
library(cowplot)
library(foreach)

snpgdsVCF2GDS(vcf.fn = "/mnt/pricey_2/priscilla/hybrid/hc/hs.hc.both.keep.linenames.vcf",  method=c("biallelic.only"), snpfirstdim = FALSE, out.fn="/mnt/pricey_2/priscilla/hybrid/hc/hs.hc.both.keep.linenames.vcf.gds" )

snpgdsVCF2GDS(vcf.fn = "/mnt/pricey_2/priscilla/hybrid/hc/hs.hc.both.keep.linenames.vcf",  method=c("biallelic.only"), snpfirstdim = FALSE, out.fn="/mnt/pricey_2/priscilla/hybrid/hc/hs.hc.both.keep.linenames.vcf.gds" )

genofile <- snpgdsOpen("/mnt/pricey_2/priscilla/hybrid/hc/hs.hc.both.keep.linenames.vcf.gds" )

# snpset <- snpgdsLDpruning(genofile,
#                           ld.threshold=0.2,
#                           slide.max.bp = 5000,
#                           autosome.only=FALSE,
#                           maf=.05)
# 
# snpset.id <-unlist(snpset)

a<-snpgdsSNPList(genofile)
info<-data.table(snp.id=a$snp.id,
                 chr=a$chromosome,
                 pos=a$pos,
                 freq=a$afreq)

pcaOut <- foreach(i=c("all", "2L", "2R", "3L" ,"3R", "X"))%do%{
    print(i)
    if(i=="all") snpset.id.use <- info$snp.id
    if(i!="all") snpset.id.use <- info[chr==i, snp.id]
    pca.temp <- snpgdsPCA(genofile, snp.id=snpset.id.use, autosome.only=FALSE, num.thread=10)
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
                         
                         sample.id=pca.temp$sample,
                         chrs=i)
    pca.dt
}
pcaOut <- rbindlist(pcaOut)

lines<-fread("/mnt/pricey_2/priscilla/all_line_data.txt")
lines[,sample.id:=strain]
pcaOut=merge(pcaOut, lines, by="sample.id")

ggplot(pcaOut[chrs=="all"], aes(x=pc1, y=pc2, color=geography))+geom_point()
ggplot(pcaOut[chrs=="all"], aes(x=pc3, y=pc2, color=geography))+geom_point()

write.table(pcaOut, "/mnt/pricey_2/priscilla/revisions/parents_PCA.txt", sep="\t", quote=F, row.names=F)

pcaOut=fread("/mnt/pricey_2/priscilla/revisions/parents_PCA.txt")


#now do spring vs fall only

seasons<-lines[grepl("12LN", strain), strain]

pcaOut <- foreach(i=c("all", "2L", "2R", "3L" ,"3R", "X"))%do%{
    print(i)
    if(i=="all") snpset.id.use <- info$snp.id
    if(i!="all") snpset.id.use <- info[chr==i, snp.id]
    pca.temp <- snpgdsPCA(genofile, snp.id=snpset.id.use, sample.id=seasons, autosome.only=FALSE, num.thread=10)
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
                         
                         sample.id=pca.temp$sample,
                         chrs=i)
    pca.dt
}
pcaOut <- rbindlist(pcaOut)

pcaOut=merge(pcaOut, lines, by="sample.id")
write.table(pcaOut, "/mnt/pricey_2/priscilla/revisions/parents_season_PCA.txt", sep="\t", quote=F, row.names=F)


ns<-lines[!grepl("DGRP", strain), strain]

pcaOut <- foreach(i=c("all", "2L", "2R", "3L" ,"3R", "X"))%do%{
    print(i)
    if(i=="all") snpset.id.use <- info$snp.id
    if(i!="all") snpset.id.use <- info[chr==i, snp.id]
    pca.temp <- snpgdsPCA(genofile, snp.id=snpset.id.use, sample.id=ns, autosome.only=FALSE, num.thread=10)
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
                         
                         sample.id=pca.temp$sample,
                         chrs=i)
    pca.dt
}
pcaOut <- rbindlist(pcaOut)

pcaOut=merge(pcaOut, lines, by="sample.id")
write.table(pcaOut, "/mnt/pricey_2/priscilla/revisions/parents_northsouth_PCA.txt", sep="\t", quote=F, row.names=F)





library(GGally)
ggpairs(pcaOut, columns=c(2:11), aes(color=geography))

pcaOut[geography=="Maine", loc:="North"]
pcaOut[geography=="Penn-spring", loc:="North"]
pcaOut[geography=="Penn-fall", loc:="North"]
pcaOut[geography=="New York", loc:="North"]
pcaOut[geography=="North Carolina", loc:="South"]
pcaOut[geography=="Southeast", loc:="South"]
pcaOut[geography=="Carribean", loc:="South"]

ggpairs(pcaOut, columns=c(2:11), aes(color=loc))
ggpairs(pcaOut[geography=="Penn-spring"|geography=="Penn-fall"], columns=c(2:11), aes(color=geography))


#moved all these files to sammas storage 1/1/2020




