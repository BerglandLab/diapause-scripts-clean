library(data.table)
library(foreach)
library(data.table)
library(foreach)

import_data<-function(fn){
    if (!file.size(fn) == 0) {
        dat<-fread(fn, header=TRUE)
        #rint(fn)
        sample=strsplit(fn,split=".haps")[[1]]
        print(sample)
        dat[, sample.id:=rep(sample, dim(dat)[1])]
        return(dat)
    }
}


#read in simulated reconstructed paths
all.sim.rec<-foreach(swarm=c("A", "B"), .combine="rbind")%do%{
    wd=paste0("/scratch/pae3g/genome-reconstruction/", swarm)
    print(wd)
    setwd(wd)
    file_list <- list.files(pattern="34F.*estimate.14.haps")
    keep <-foreach(fn=file_list, .combine="rbind") %do% {
        import_data(fn)
    }
    return(keep)
}

all.sim.rec[, swarm:=tstrsplit(sample.id, split="[.]")[[4]]]
all.sim.rec[, ind:=tstrsplit(sample.id, split="[.]")[[3]]]
all.sim.rec[, gen:=substring(tstrsplit(sample.id, split="_")[[2]], 1,1)]
all.sim.rec.melt=melt(all.sim.rec, id.vars = c("chromosome", "start", "stop", "sample.id", "swarm", "gen"), measure.vars=c("par1", "par2"), variable.name = "haplotype", value.name="line")

all.sim.rec.melt[, path.length:=stop-start]
all.sim.rec.melt=all.sim.rec.melt[order(sample.id, chromosome, haplotype,start)]
all.sim.rec.melt[,run:=rleid(sample.id, chromosome, haplotype, line, gen)]


write.table(all.sim.rec.melt, "/scratch/pae3g/genome-reconstruction/all_simulated_reconstructed_haplotypes_melted.txt", sep="\t", quote=F, row.names=F)