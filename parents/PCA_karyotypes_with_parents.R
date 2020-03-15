
#on rivanna

library(gdsfmt)
library(SNPRelate)
library(data.table)
library(cowplot)
library(foreach)
library(lattice)
library(tidyr)
library(SeqArray)
library(stringr)
library(doMC)
registerDoMC(20)

#first, get a set of informative snps from parents
#make gds of full parental genotype file 

#snpgdsVCF2GDS("/mnt/pricey_2/priscilla/hybrid/hc/hs.hc.snp.99.9.both.linenames.vcf",  "/mnt/sammas_storage/bergland-lab/Priscilla/hs.hc.snp.99.9.both.linenames.vcf.gds", method="biallelic.only")

#make a new gds of final hybrid swarm
#snpgdsVCF2GDS("/scratch/pae3g/oldscratch_recovered/genome-reconstruction/final2.vcf", "/scratch/pae3g/revisions/final2.remade.vcf.gds", method="biallelic.only")

#try merging
#snpgdsCombineGeno(gds.fn = c("/nv/vol186/bergland-lab/Priscilla/hs.hc.snp.99.9.both.linenames.vcf.gds", "/scratch/pae3g/revisions/final2.remade.vcf.gds"), out.fn = "/scratch/pae3g/revisions/hybrid_parent_allsnp_merge.vcf.gds", same.strand = T)


genofile <- snpgdsOpen("/scratch/pae3g/revisions/hybrid_parent_allsnp_merge.vcf.gds" , allow.fork=T)

#get informative
snpset <- snpgdsLDpruning(genofile,
                          ld.threshold=0.2,
                          slide.max.bp = 5000,
                          autosome.only=FALSE,
                          maf=.05)

snpset.id <-unlist(snpset)

a<-snpgdsSNPList(genofile)
info<-data.table(snp.id=a$snp.id,
                 chr=a$chromosome,
                 pos=a$pos,
                 freq=a$afreq)

snpset <- snpgdsLDpruning(geno,
                          ld.threshold=0.2,
                          slide.max.bp = 5000,
                          autosome.only=FALSE,
                          maf=.05)

snpset.id <-unlist(snpset)

 #empty list to save PCs
 pcaVars=list()

 #calculate PCs
 pcaOut <- foreach(i=c("all", "chr2L", "chr2R", "chr3L" ,"chr3R", "chrX"))%do%{
     print(i)
     if(i=="all") snpset.id.use <- snpset.id
     if(i!="all") snpset.id.use <- unlist(snpset[which(names(snpset)==i)])
     print(length(snpset.id.use))
     pca.temp <- snpgdsPCA(genofile ,snp.id=snpset.id.use, autosome.only=FALSE, num.thread=10)
     pcaVars[[i]]=pca.temp$varprop
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
 
 metadata=fread("/scratch/pae3g/oldscratch_recovered/final_reconstruction/PAE_AOB_library_metadata_all.txt", header=FALSE)
 
 setnames(metadata, c("sample.id", "generation", "library", "actual.cage"))
 
 pcaOut2=merge(pcaOut, metadata, by="sample.id", all.x=T)

 
 lines<-fread("/nv/vol186/bergland-lab/Priscilla/all_line_data.txt")
 setnames(lines, "strain", "sample.id")
 pcaOut2<-merge(pcaOut2, lines[, .(sample.id, geography, population)], all.x=T)
 
 pcaOut2[, type:=ifelse(grepl("PAE", sample.id), "hybrid", "parent")]
 pcaOut2[,group:=ifelse(type=="hybrid", paste(type, actual.cage, sep="-"),paste(type, population, sep="-"))]

save(pcaOut2, file="/scratch/pae3g/revisions/hybrid_parent_PCA.Rdata")

#get karyotypes 


 snp.dt <- data.table(snp.id=read.gdsn(index.gdsn(genofile, "snp.id")),
                      chr=read.gdsn(index.gdsn(genofile, "snp.chromosome")),
                      pos=read.gdsn(index.gdsn(genofile, "snp.position")),
                      alleles=read.gdsn(index.gdsn(genofile, "snp.allele")))

 ### load inversion markers
 fixedInv <- read.delim("/scratch/pae3g/oldscratch_recovered/evolution/inv_fixed_coord.txt", header=F, sep=" ", as.is=T)
 fixedInv <- as.data.frame(t(fixedInv))
 fixedInv$chr <- substr(fixedInv[,1], 4, 5)
 fixedInv$pos <- as.numeric(as.character(fixedInv[,2]))

 fixedInv <- as.data.table(fixedInv)

 ### merge
 setkey(snp.dt, chr, pos)
 setkey(fixedInv, chr, pos)
 snp.inv <- merge(snp.dt, fixedInv, all.y=T)

 ### get allelic status
 genomat.orig <- snpgdsGetGeno(genofile, snp.id=na.omit(snp.inv$snp.id), with.id=T)
 genomat <- data.table(genotype=expand.grid(as.matrix(genomat.orig$genotype))$Var1,
                       sample.id=rep(genomat.orig$sample.id, length(genomat.orig$snp.id)),
                       snp.id=rep(genomat.orig$snp.id, each=length(genomat.orig$sample.id)))

 setkey(genomat, snp.id)


 ### merge back with inversion names
 setkey(snp.inv, snp.id)
 genomat <- merge(genomat, snp.inv)

 ### split alt and ref
 genomat[,ref:=do.call("rbind", strsplit(alleles, "/"))[,1]]
 genomat[,alt:=do.call("rbind", strsplit(alleles, "/"))[,2]]

 ### throw out remaining indels
 genomat <- genomat[nchar(alt)%in%c(1,3)]

 ### is 'inversion' allele in alt?
 genomat[nchar(alt)==3, invCheck := apply(unlist(sapply(V3, grepl, alt)), 2, any)]
 genomat[nchar(alt)==1, invCheck := V3==alt]

 genomat <- genomat[invCheck==T]

 ### deal with monomorphic sites from fixedInv
 setkey(genomat, chr, pos)
 genomat.sites <- genomat[!duplicated(genomat)][,c("chr", "pos"), with=F]
 genomat.sites[,poly:=T]

 setkey(genomat.sites, chr, pos)
 setkey(fixedInv, chr, pos)

 gf <- merge(genomat.sites, fixedInv, all.x=T, all.y=T)
 gf[is.na(poly), poly:=F]

 #### okay, so ~150 sites from Martin's data are monomorphic in our data
 #### let's assume that they are all reference alleles

 mono <- data.table(snp.id=NA,
                    genotype=2,
                    sample.id=rep(unique(genomat$sample.id), each=sum(!gf$poly)),
                    chr=rep(gf[poly==F]$chr, length(unique(genomat$sample.id))),
                    pos=rep(gf[poly==F]$pos, length(unique(genomat$sample.id))),
                    alleles=NA,
                    V1=rep(gf[poly==F]$V1, length(unique(genomat$sample.id))),
                    V2=rep(gf[poly==F]$V2, length(unique(genomat$sample.id))),
                    V3=rep(gf[poly==F]$V3, length(unique(genomat$sample.id))),
                    ref=NA,
                    alt=NA,
                    invCheck=T)



 genomat <- rbind(genomat, mono)


 ### make karyotype calls
 ### tabulate calls
 genomat.ag <- genomat[,list(n.std=sum(genotype==2, na.rm=T),
                             n.het=sum(genotype==1, na.rm=T),
                             n.inv=sum(genotype==0, na.rm=T)),
                       list(sample.id, inversion=V1)]
 genomat.ag <- as.data.table(gather(genomat.ag, kary, n, n.std:n.inv))

 ### make most-likely call
 kary <- genomat.ag[,list(kary=kary[which.max(n)],
                          evidence=max(n)/sum(n),
                          n=sum(n),
                          prStd=(2*n[kary=="n.std"] + n[kary=="n.het"])/(2*sum(n))),
                    list(sample.id, inversion)]

 ### a plot
 ggplot(data=kary, aes(x=inversion, y=prStd, color=inversion)) +
     geom_boxplot() + geom_hline(yintercept=.9)

 ### are there any strains where the pr(Std) is < .99 for more than two inversions per chromosome?
 kary[,chr:=gsub("In\\(([2,3,X,L,R]{1,})\\)[A-Za-z]{1,}", "\\1", inversion)]

 kary.ag <- kary[,list(n=sum(prStd<.9)), list(chr, sample.id)]
 table(kary.ag$n)

 ### make final karyotype calls
 karyCall <- function(prStds, inv, cutVal=.75) {
     if(all(prStds>=cutVal)) {
         o <- "std"
     } else if (sum(prStds <= (1-cutVal))==1) {
         o <- as.character(inv[prStds<=(1-cutVal)])
     } else if (sum(prStds>=(1-cutVal) & prStds<=cutVal)==2) {
         o <- paste(inv[prStds>=(1-cutVal) & prStds<=cutVal], collapse=";")
     } else {
         o <- paste("std", as.character(inv[prStds>=(1-cutVal) & prStds<=cutVal]), sep=";")
     }

     return(o)
 }

 kary.ag <- kary[!is.nan(prStd), list(chr2L=karyCall(prStds=prStd[chr=="2L"], inv=inversion[chr=="2L"]),
                                      chr2R=karyCall(prStds=prStd[chr=="2R"], inv=inversion[chr=="2R"]),
                                      chr3L=karyCall(prStds=prStd[chr=="3L"], inv=inversion[chr=="3L"]),
                                      chr3R=karyCall(prStds=prStd[chr=="3R"], inv=inversion[chr=="3R"])),
                 list(sample.id)]


save(kary.ag, file="/scratch/pae3g/revisions/parents_hybrids_karyotype_calls.Rdata") #note that this file seems to miss some of the parental inversions...


b<-merge(kary.ag)