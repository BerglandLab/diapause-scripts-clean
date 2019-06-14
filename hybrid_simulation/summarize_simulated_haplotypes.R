
library(data.table)
library(foreach)
library(cowplot)

#pull in the population files for 5 individuals from 100 populations of each swarm/generation

a<-foreach(swarm=c("A", "B"), .combine="rbind") %do%{
    b<-foreach(gen=c(4:5), .combine="rbind")%do%{
        c<-foreach(rep=c(1:100), .combine="rbind")%do%{
            d<-fread(paste0("/nv/vol186/bergland-lab/Priscilla/genome-reconstruction/output/cage",swarm, "/34F_", gen, "G_0.5X_10000N_rep", rep, "/population.haps"))
            d<-d[ind<=5]
            d[,swarm:=swarm]
            d[,gen:=gen]
            d[,rep:=rep]
            return(d)
            
        }
        return(c)
    }
    return(b)
}

#a has simulated haplotypes for 5000 individuals for each combination of swarm and generation
#make ID column
a[,id:=paste(swarm, gen, rep, ind, sep="_")]

write.table(a, "/scratch/pae3g/genome-reconstruction/simulated_actual_haplotyp> .txt", sep="\t", quote=F, row.names=F)