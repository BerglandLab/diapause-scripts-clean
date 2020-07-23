library(data.table)

keep.paths<-fread("/scratch/pae3g/evolution/all_reconstructed_haplotypes_melted.txt")

keep.paths[path.length<1000000, line:="UNKNOWN"]
paths.sum<-keep.paths[,.(prop.unknown=sum(path.length[line=="UNKNOWN"])/sum(path.length), n.chr.imp=length(unique(chromosome))), .(sample.id, swarm)]

phenos<-fread("/nv/vol186/bergland-lab/Priscilla/phenos_wolbachia.txt")

phenos<-merge(paths.sum[,.(sample.id, prop.unknown, n.chr.imp)], phenos, by="sample.id")

write.table(phenos, "/scratch/pae3g/revisions/phenos_wolbachia_missinggeno.txt", row.names=F, quote=F, sep="\t")