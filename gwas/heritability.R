#make covar file in R
library(data.table)
library(foreach)

phenos <- fread("/scratch/pae3g/oldscratch_recovered/phenos/phenos_062018.txt")
 setkey(phenos, id)


 files<-fread("/scratch/pae3g/genome-reconstruction/new_permutation_input5.txt")
 
 #read in files
foreach( perm=c(files$V5[1:1000], 0))%do%{
     ### generate permuted phenotype, qcov, and covfiles for GATC

     print(perm)
     seed<-files[V5==perm,V6]
     print(seed)
     set.seed(seed)
     phenos[swarm=="A", id:=sample(id)]
     phenos[swarm=="B", id:=sample(id)]
     write.table(phenos[,.(id, id, diapause.bin9, diapause.bin, Eggs, eggP)], paste0("/scratch/pae3g/revisions/gcta/plink_phenos_perm", perm, ".txt"), quote=F, sep="\t", row.names=F, col.names=F)
     
 
#run heritabiltity on workshtation with bash script

#read in results

library(data.table)
library(foreach)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(ggbeeswarm)


y<-foreach(perm=c(0,1000:1999), .errorhandling="remove")%do%{
    h<-fread(paste0("/mnt/pricey_2/priscilla/gcta/gcta_", perm, ".hsq"), fill=T)
    h[,perm:=perm]
    return(h)
}

y<-rbindlist(y)

write.table(y,"/mnt/pricey_2/priscilla/gcta_hsq_1000perms.txt", quote=F, sep="\t", row.names=F)

h.plot<-ggplot(y[perm!=0&(Source=="V(G)"|Source=="V(e)"|Source=="V(G)/Vp")], aes(x=Source, y=Variance))+
    geom_quasirandom(color="grey", varwidth = TRUE, size=.5, method = "smiley")+
    geom_point(data=y[perm==0&(Source=="V(G)"|Source=="V(e)"|Source=="V(G)/Vp")], aes(x=Source, y=Variance), color="lightseagreen", size=2)+
    labs(y=NULL, x=NULL, title="Stage 9")+
    scale_x_discrete(labels=c(expression("V"["e"]), expression("V"["g"]), expression("h"^2)))+
    lims(y=c(0, .18))+
    geom_errorbar(data=y[perm==0&(Source=="V(G)"|Source=="V(e)"|Source=="V(G)/Vp")], aes(x=Source, ymin=Variance-1.96*SE, ymax=Variance+1.96*SE), color="lightseagreen", width=0.05)






z<-foreach(perm=c(0,1000:1999), .errorhandling="remove")%do%{
    h<-fread(paste0("/mnt/pricey_2/priscilla/gcta/gcta_st7", perm, ".hsq"), fill=T)
    h[,perm:=perm]
    return(h)
}

z<-rbindlist(z)

write.table(z,"/mnt/pricey_2/priscilla/gcta_hsq_1000perms_st7.txt", quote=F, sep="\t", row.names=F)

