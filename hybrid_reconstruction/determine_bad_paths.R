library(data.table)
library(cowplot)
library(foreach)

#some functions

#imports haplotype path files and pulls individual id out of file name
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

#function that melts diploid haplotypes into single haplotypes and calculates start and stop postiion of each haplotype
meltPaths=function(obj){
    #read in metadata for all samples
    metadata=fread("/scratch/pae3g/final_reconstruction/PAE_AOB_library_metadata_all.txt", header=FALSE)
    setnames(metadata, c("sample.id", "generation", "library", "swarm"))
    #merge in metadata to path object
    all.paths=merge(obj, metadata, by="sample.id")
    #get number-> line translation
    #parent_line_by_swarm.txt translates the "p###" reconstruction notation to actual line names
    swarms=fread("/scratch/pae3g/final_reconstruction/parent_line_by_swarm.txt", header=FALSE)
    names(swarms)=c("line", "num.68", "numeric", "swarm")
    #melt file to long format so there is a single haplotype per line
    all.paths.melt=melt(all.paths, id.vars = c("chromosome", "start", "stop", "sample.id", "swarm"), measure.vars=c("par1", "par2"), variable.name = "haplotype", value.name="numeric")
    #calculate path length
    all.paths.melt[, path.length:=stop-start]
    #merge in line info
    all.paths.melt=merge(all.paths.melt, swarms, by=c("swarm", "numeric"))
    #get rid of extra column
    all.paths.melt[,num.68:=NULL]
    #order by id, chr, position
    all.paths.melt=all.paths.melt[order(sample.id, chromosome, haplotype,start)]
    #assign a run id for all lines with the same sample number, chromosome, haplotype, and founder line
    all.paths.melt[,run:=rleid(sample.id, chromosome, haplotype, line)]
    #find maximum and minimum (start and stop) for consecutive paths
    all.paths.cons=all.paths.melt[,.(cons.start=min(start), cons.stop=max(stop)), .(sample.id, chromosome, haplotype, line, swarm, run)]
    #calculate new path lengths
    all.paths.cons[,path.length:=cons.stop-cons.start]
    #merge in metadata
    all.paths.cons=merge(all.paths.cons, metadata, by=c("sample.id", "swarm"))
    return(all.paths.cons[order(sample.id, chromosome, haplotype, cons.start)])
}

#pull in all data for reconstruction run with 14 possible founders
keep.14 <- foreach(swarm=c("A", "B"), .combine="rbind") %do%{
    setwd(paste("/scratch/pae3g/genome-reconstruction/", swarm, sep=""))
    #get all estimated haplotype file names
    file_list <- list.files(pattern="PAE.*estimate.14.haps")
    #import all data
    swarm.list<-foreach(fn=file_list, .combine="rbind") %do% {
        import_data(fn)
    }
    return(swarm.list)
}

#melt into single haplotypes
all.paths.14<-meltPaths(keep.14)

#save path file
write.csv(all.paths.14, "/scratch/pae3g/final_reconstruction2/all_paths.csv")

#calculate the proportion of paths that are short (<1 Mb)
sum(as.numeric(all.paths.14[path.length<1000000, path.length]))/sum(as.numeric(all.paths.14[,path.length]))#0.89%

#look at number of recombinations per chromosome on collapsed data
r.14=all.paths.14[,.(Nrecombs.new.14=.N-1), .(chromosome, sample.id)]

#deal with missing data to find paths to mask
#assign all paths <1 Mb to "unknown"
all.paths.14[path.length<1000000, line:="UNKNOWN"]


#calculate some summaries

#caclulate the proportion of missing genotypes for each individual
mis=all.paths.14[,.(prop.missing=sum(path.length[line=="UNKNOWN"])/sum(path.length)), .(sample.id)]

#caculate proportion missing genotypes for each individual chromosome arm
mis.chr=all.paths.14[,.(prop.missing=sum(path.length[line=="UNKNOWN"])/sum(path.length)), .(sample.id, chromosome)]

#merge in missing path data

all.paths.14=merge(all.paths.14, mis.chr, by=c("sample.id", "chromosome"))

#if a chromosme is >20% missing data, assign the entire chromosome as missing/unknown
all.paths.14[prop.missing>.2, line:="UNKNOWN"]

#caculate new percentage missing data
sum(as.numeric(all.paths.14[line=="UNKNOWN", path.length]))/sum(as.numeric(all.paths.14[,path.length]))#1.45% now

#pull out all unknown/missing data
bad.paths=all.paths.14[line=="UNKNOWN"]

#make file to be used in vcf masking script
write.csv(bad.paths, "/scratch/pae3g/final_reconstruction2/bad_paths.csv")


#now use this bad path file to recreate vcf/gds. any short paths (under 1000000 bp) are unknown, and any chromosomes with more than 20% short is also blanked out.


#continue processing to caculate number of recombinants -- bridge over missing data <1 Mb that has same haplotype on either side:

#drop all paths that are still < 1 Mb even after merging 
real.merge<-all.paths.14[path.length>1000000]

#reassign new rleids
real.merge[,run:=rleid(sample.id, chromosome, haplotype, line)]

#calculate new starts and stops
real.merge2=real.merge[,.(cons.start=min(cons.start), cons.stop=max(cons.stop)), .(sample.id, chromosome, haplotype, line, swarm, run,  gen)]

#calculate new lengths
real.merge2[,path.length:=cons.stop-cons.start]

#count recombinations in each dataset (# of recomb = # haplotypes -1 )
real.merge2.r<-real.merge2[,.(R=.N-1),.(sample.id,chromosome, swarm, gen)]
