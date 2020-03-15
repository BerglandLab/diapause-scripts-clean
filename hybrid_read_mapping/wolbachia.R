library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

files<-fread("/nv/vol186/bergland-lab/Priscilla/hybrid_bam/files.txt", header=F)$V1

b<-foreach (fn=files, .errorhandling="remove")%do%{
    a<-fread(paste0("/nv/vol186/bergland-lab/Priscilla/hybrid_bam/", fn, ".stats"), header=F)
    a[,file:=fn]
    return(a)
}

b<-rbindlist(b)

setnames(b, c("V1", "V2", "V3", "V4"),c("chr", "length", "mapped", "unmapped" ))

b.sum<-b[,.(chrom.reads=sum(mapped[chr%in%c("2L", "2R", "3L", "3R", "X")]),
            wol.reads=mapped[chr=="gi|42519920|ref|NC_002978.6|"]), .(file)]

b.sum[,prop.wolb:=log10(wol.reads/(wol.reads+chrom.reads))]

ggplot(b.sum, aes(x=prop.wolb))+geom_histogram()+labs(x="log10(wolbachia reads: chromosomal reads)")

b.sum[,id:=tstrsplit(file, split=".merged")[[1]]]

phenos <- fread("/scratch/pae3g/oldscratch_recovered/phenos/phenos_062018.txt")

b.sum<-merge(b.sum, phenos, by="id")

b.sum[,wolbachia:=ifelse(prop.wolb<(-3), F, T)]

table(b.sum$diapause.bin9,b.sum$wolbachia)
table(b.sum$diapause.bin,b.sum$wolbachia)

fisher.test(table(b.sum$diapause.bin9,b.sum$wolbachia))
fisher.test(table(b.sum$diapause.bin,b.sum$wolbachia))


summary(glm(diapause.bin9~temp.rack.cal+photoperiod+swarm+generation+wolbachia, data=b.sum, family="binomial"))

summary(glm(diapause.bin~temp.rack.cal+photoperiod+swarm+generation+wolbachia, data=b.sum, family="binomial"))


write.table(b.sum, "/nv/vol186/bergland-lab/Priscilla/phenos_wolbachia.txt", sep="\t", quote=F, row.names=F)


#look at parents

files<-fread("/nv/vol186/bergland-lab/Priscilla/hs_parents_mapped/files.txt", header=F)$V1

b<-foreach (fn=files, .errorhandling="remove")%do%{
    a<-fread(paste0("/nv/vol186/bergland-lab/Priscilla/hs_parents_mapped/", fn, ".stats"), header=F)
    a[,file:=fn]
    return(a)
}

b<-rbindlist(b)

setnames(b, c("V1", "V2", "V3", "V4"),c("chr", "length", "mapped", "unmapped" ))

b.sum<-b[,.(chrom.reads=sum(mapped[chr%in%c("2L", "2R", "3L", "3R", "X")]),
            wol.reads=mapped[chr=="gi|42519920|ref|NC_002978.6|"]), .(file)]

b.sum[,prop.wolb:=log10(wol.reads/(wol.reads+chrom.reads))]

ggplot(b.sum, aes(x=prop.wolb))+geom_histogram()+labs(x="log10(wolbachia reads: chromosomal reads)")

b.sum[,line:=tstrsplit(file, split="[.]")[[1]]]

b.sum[,wolbachia:=ifelse(prop.wolb<(-3), F, T)]

write.csv(b.sum, file="/scratch/pae3g/revisions/hybrid_parents_wolbachia.csv")
