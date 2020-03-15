library(SeqArray)
library(data.table)

seqVCF2GDS(vcf.fn="/mnt/pricey_2/priscilla/hs.hc.snp.99.9.both.linenames.ann.vcf", out.fn="/mnt/pricey_2/priscilla/hs.hc.snp.99.9.both.linenames.ann.vcf.seq.gds")

genofile <- seqOpen("/mnt/pricey_2/priscilla/hs.hc.snp.99.9.both.linenames.ann.vcf.seq.gds")


snp.dt <- data.table(variant.id=seqGetData(genofile, "variant.id"))
tmp <- seqGetData(genofile, "annotation/info/ANN")
len1 <- tmp$length
len2 <- tmp$data

snp.dt1 <- data.table(len=rep(len1, times=len1),
                      ann=len2,
                      id=rep(snp.dt$variant.id, times=len1))

snp.dt1[,class:=tstrsplit(snp.dt1$ann,"\\|")[[2]]]
snp.dt1[,gene:=tstrsplit(snp.dt1$ann,"\\|")[[4]]]


# Collapsing additional annotations to original SNP vector length
snp.dt1.an <- snp.dt1[,list(n=length(class), col.class= paste(class, collapse=","), col.gen=paste(gene, collapse=",")),
                      list(variant.id=id)]

snp.dt1.an[,col.class:=tstrsplit(snp.dt1.an$col.class,"\\,")[[1]]]
snp.dt1.an[,col.gene:=tstrsplit(snp.dt1.an$col.gene,"\\,")[[1]]]

# merge

m <- merge(snp.dt, snp.dt1.an, by="variant.id")

m[,chr:=seqGetData(genofile, "chromosome")]
m[,pos:=seqGetData(genofile, "position")]


save(m, file="/mnt/pricey_2/priscilla/snp_annotations.Rdat")
# merge


