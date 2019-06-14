library(data.table)
library(foreach)
library(cowplot)
library(ggbeeswarm)

setwd("/scratch/pae3g/genome-reconstruction/accuracy")

fn<-list.files(pattern=".dat")

acc<-foreach(file=fn, .combine="rbind", .errorhandling = "remove")%do%{
    return(fread(file))
}

acc[,swarm:=tstrsplit(ind_id, split="[.]")[[4]]]

acc[,group:=paste("Pop.", swarm, "gen.", nGenerations, sep=" ")]

write.table(acc, "/scratch/pae3g/genome-reconstruction/accuracy.txt", quote=F, row.names=F, sep="\t" )