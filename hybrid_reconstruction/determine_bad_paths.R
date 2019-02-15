library(data.table)
library(cowplot)
library(foreach)

#some functions

import_data<-function(fn){
  if (!file.size(fn) == 0) {
    dat<-fread(fn, header=TRUE)
    #rint(fn)
    sample=strsplit(fn,split="[.]")[[1]][1]
    print(sample)
    dat[, sample.id:=rep(sample, dim(dat)[1])]
    return(dat)
  }
}

meltPaths=function(obj){
  metadata=fread("/scratch/pae3g/final_reconstruction/PAE_AOB_library_metadata_all.txt", header=FALSE)
  setnames(metadata, c("sample.id", "generation", "library", "swarm"))
  #merge in metadata
  all.paths=merge(obj, metadata, by="sample.id")
  #get number-> line translation
  swarms=fread("/scratch/pae3g/final_reconstruction/parent_line_by_swarm.txt", header=FALSE)
  names(swarms)=c("line", "num.68", "numeric", "swarm")
  #melt file so single haplotype per line
  all.paths.melt=melt(all.paths, id.vars = c("chromosome", "start", "stop", "sample.id", "swarm"), measure.vars=c("par1", "par2"), variable.name = "haplotype", value.name="numeric")
  all.paths.melt[, path.length:=stop-start]
  all.paths.melt=merge(all.paths.melt, swarms, by=c("swarm", "numeric"))
  all.paths.melt[,num.68:=NULL]
  all.paths.melt=all.paths.melt[order(sample.id, chromosome, haplotype,start)]
  all.paths.melt[,run:=rleid(sample.id, chromosome, haplotype, line)]
  all.paths.cons=all.paths.melt[,.(cons.start=min(start), cons.stop=max(stop)), .(sample.id, chromosome, haplotype, line, swarm, run)]
  all.paths.cons[,path.length:=cons.stop-cons.start]
  all.paths.cons=merge(all.paths.cons, metadata, by=c("sample.id", "swarm"))
  return(all.paths.cons[order(sample.id, chromosome, haplotype, cons.start)])
}


keep.14 <- foreach(swarm=c("A", "B"), .combine="rbind") %do%{
  setwd(paste("/scratch/pae3g/genome-reconstruction/", swarm, sep=""))
  file_list <- list.files(pattern="PAE.*estimate.14.haps")
  swarm.list<-foreach(fn=file_list, .combine="rbind") %do% {
    import_data(fn)
  }
  return(swarm.list)
}

#melt into single haplotypes
all.paths.14<-meltPaths(keep.14)

#all 2823 samples worked


#plot haplotype sizes
ggplot()+geom_density(data=all.paths.14, aes(x=log10(path.length), color=as.factor(generation)))

#calculate the proportion of paths that are short
sum(as.numeric(all.paths.14[path.length<1000000, path.length]))/sum(as.numeric(all.paths.14[,path.length]))#0.89%



#look at number of recombinations per chromosome
r.14=all.paths.14[,.(Nrecombs.new.14=.N-1), .(chromosome, sample.id)]


#deal with missing data to find paths to mask
all.paths.14[path.length<1000000, line:="UNKNOWN"]
write.csv(all.paths.14, "/scratch/pae3g/final_reconstruction2/all_paths.csv")


mis=all.paths.14[,.(prop.missing=sum(path.length[line=="UNKNOWN"])/sum(path.length)), .(sample.id)]
mis.chr=all.paths.14[,.(prop.missing=sum(path.length[line=="UNKNOWN"])/sum(path.length)), .(sample.id, chromosome)]

failed.chr=mis.chr[prop.missing>.2] #122 failed chromosomes

all.paths.14=merge(all.paths.14, mis.chr, by=c("sample.id", "chromosome"))
all.paths.14[prop.missing>.2, line:="UNKNOWN"]

sum(as.numeric(all.paths.14[line=="UNKNOWN", path.length]))/sum(as.numeric(all.paths.14[,path.length]))#1.45% now
bad.paths=all.paths.14[line=="UNKNOWN"]

write.csv(bad.paths, "/scratch/pae3g/final_reconstruction2/bad_paths.csv")


#now use this bad path file to recreate vcf/gds. any short paths (under 1000000 bp) are unknown, and any chromosomes with more than 20% short is also blanked out.
