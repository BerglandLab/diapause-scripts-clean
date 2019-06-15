
#look at missing rate and heterozygosity in each founding population to compare to missing data.

#first need gds files of each populations

library(SNPRelate)
library(foreach)
library(data.table)

A <- snpgdsOpen("/mnt/pricey_2/priscilla/hybrid/hc/hs.hc.A.vcf.gds")

snp.dt <- data.table(snp.id=read.gdsn(index.gdsn(A, "snp.id")),
                     chr=read.gdsn(index.gdsn(A, "snp.chromosome")),
                     pos=read.gdsn(index.gdsn(A, "snp.position")),
                     alleles=read.gdsn(index.gdsn(A, "snp.allele")))


#get missing data

A.mis<-as.data.table(do.call(cbind,snpgdsSNPRateFreq(A, with.snp.id=TRUE)))

#make table of rate of heterozygosity
par.geno=snpgdsGetGeno(A, with.id=TRUE)
genos=par.geno$genotype
rownames(genos)=par.geno$sample.id
colnames(genos)=par.geno$snp.id
genos2=t(genos)
genos3=as.data.table(genos2)
genos3[,snp.id:=par.geno$snp.id]
genos4=melt(genos3,id.var="snp.id", measure.vars=names(genos3)[1:34])
A.het=genos4[,.(num.het=sum(value==1, na.rm=TRUE), num.ref=sum(value==2, na.rm=TRUE), num.alt=sum(value==0, na.rm=TRUE), num.missing=sum(is.na(value))), .(snp.id)]



A.het[,prop.het:=num.het/34]
hist(A.het$prop.het, breaks=34)
#total number of snps to drop =578718

A.het=merge(A.het, snp.dt, by="snp.id")
A.het[,start:=pos-1]

A.keep=A.het[(num.het+num.missing)<3 & num.alt>0 & num.ref>0] #467,104

A.keep=merge(A.keep, snp.dt, by="snp.id")
A.keep[,start:=pos-1]
A.keep.bed=A.keep[,.(chr, start, pos)]

write.table(A.keep.bed, "/mnt/pricey_2/priscilla/hybrid/hc/hs.hc.A.keep.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)


#same thing for B

B <- snpgdsOpen("/mnt/pricey_2/priscilla/hybrid/hc/hs.hc.B.vcf.gds")

snp.dt <- data.table(snp.id=read.gdsn(index.gdsn(B, "snp.id")),
                     chr=read.gdsn(index.gdsn(B, "snp.chromosome")),
                     pos=read.gdsn(index.gdsn(B, "snp.position")),
                     alleles=read.gdsn(index.gdsn(B, "snp.allele")))

#get missing data

B.mis<-as.data.table(do.call(cbind,snpgdsSNPRateFreq(B, with.snp.id=TRUE)))

#make table of rate of heterozygosity
par.geno=snpgdsGetGeno(B, with.id=TRUE)
genos=par.geno$genotype
rownames(genos)=par.geno$sample.id
colnames(genos)=par.geno$snp.id
genos2=t(genos)
genos3=as.data.table(genos2)
genos3[,snp.id:=par.geno$snp.id]
genos4=melt(genos3,id.var="snp.id", measure.vars=names(genos3)[1:34])
B.het=genos4[,.(num.het=sum(value==1, na.rm=TRUE), num.ref=sum(value==2, na.rm=TRUE), num.alt=sum(value==0, na.rm=TRUE), num.missing=sum(is.na(value))), .(snp.id)]
B.het[,prop.het:=num.het/34]
hist(B.het$prop.het, breaks=34)

sum(B.het$drop)
#total number of snps to drop =653,064
B.het=merge(B.het, snp.dt, by="snp.id")
B.het[,start:=pos-1]



B.keep=B.het[(num.het+num.missing)<3 & num.alt>0 & num.ref>0] #388,492


B.keep=merge(B.keep, snp.dt, by="snp.id")
B.keep[,start:=pos-1]
B.keep.bed=B.keep[, .(chr, start, pos)]

write.table(B.keep.bed, "/mnt/pricey_2/priscilla/hybrid/hc/hs.hc.B.keep.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)


snpgdsClose(B)