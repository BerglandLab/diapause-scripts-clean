#make covar file in R
library(data.table)
library(foreach)
library(doMC, quietly=TRUE)
registerDoMC(10)

 files<-fread("/scratch/pae3g/genome-reconstruction/universal_input.txt")
 files<-files[((V4==0&V3==1)|(V4!=0))&V6=="nonloco"&V2=="diapause.bin"]
 
 #read in files
foreach( pop=files$V1, phenotype=files$V2, perm=files$V4, seed=files$V5, draw=files$V3)%dopar%{
     ###treat phenotypes identically to mapping
    phenos <- fread("/scratch/pae3g/revisions/phenos_wolbachia_missinggeno.txt")
    setkey(phenos, id)
    filestem<- paste0("/scratch/pae3g/genome-reconstruction/final2_draw", draw, "_replaced")
    
    ### load appropriate GRM(s)	
    if (pop=="both"){
        load(paste(filestem, "_nonLOCO_Eigenstrat_allgwasSNPs.Rdat", sep=""))
    }else if (pop!="both") {
        load(paste(filestem, "_nonLOCO_Eigenstrat_allgwasSNPs_", pop, ".Rdat", sep=""))
    }
    
    #drop samples with reconstruction issues
    
    phenos <- phenos[n.chr.imp==5&prop.unknown<=0.05]
    
    #for single swarm analysis, subset phenos to only that swarm
    if(pop!="both"){
        phenos<-phenos[swarm==pop]
    }
    
    #subset sample GRM to only the samples that will be used

    #if a permutation, randomly scramble sample ids within a swarm. if perm is 0, ids will stay the same
    if(perm!=0){
        set.seed(seed)
        phenos[swarm=="A", id:=sample(id)]
        phenos[swarm=="B", id:=sample(id)]
    }

     write.table(phenos[,.(id, id, diapause.bin9, diapause.bin)], paste0("/scratch/pae3g/revisions/gcta/plink_phenos_perm", perm, "_pop", pop, ".dropmissing.txt"), quote=F, sep="\t", row.names=F, col.names=F)
     if(pop=="both"){
         write.table(phenos[,.(id, id, swarm, generation)], paste0("/scratch/pae3g/revisions/gcta/cov_gcta_perm", perm, "_pop", pop, ".dropmissing.txt"), quote=F, sep="\t", row.names=F, col.names=F)
     }else{ 
         write.table(phenos[,.(id, id, generation)], paste0("/scratch/pae3g/revisions/gcta/cov_gcta_perm", perm, "_pop", pop, ".dropmissing.txt"), quote=F, sep="\t", row.names=F, col.names=F)}
     
     write.table(phenos[,.(id, id, temp.rack.cal, photoperiod)], paste0("/scratch/pae3g/revisions/gcta/qcov_gcta_perm", perm, "_pop", pop,".dropmissing.txt"), quote=F, sep="\t", row.names=F, col.names=F)
}

     
 
     
     
# #run heritabiltity on workshtation with bash script
# 
# read in results
# 
# library(data.table)
# library(foreach)
# library(ggplot2)
# library(cowplot)
# theme_set(theme_cowplot())
# library(ggbeeswarm)
# 
# y<-list.files(path="/mnt/pricey_2/priscilla/gcta/", pattern="dropmissing.hsq")
# 
# y<-foreach(f=files, .errorhandling="remove")%do%{
#     h<-fread(f, fill=T)
#     h[,perm:=tstrsplit(f, split="_")[[3]]]
#     h[,pop:=tstrsplit(f, split="_")]
#     return(h)
# }
# 
# y<-rbindlist(y)
# 
# write.table(y,"/mnt/pricey_2/priscilla/gcta_hsq_1000perms.txt", quote=F, sep="\t", row.names=F)
# 
# h.plot<-ggplot(y[perm!=0&(Source=="V(G)"|Source=="V(e)"|Source=="V(G)/Vp")], aes(x=Source, y=Variance))+
#     geom_quasirandom(color="grey", varwidth = TRUE, size=.5, method = "smiley")+
#     geom_point(data=y[perm==0&(Source=="V(G)"|Source=="V(e)"|Source=="V(G)/Vp")], aes(x=Source, y=Variance), color="lightseagreen", size=2)+
#     labs(y=NULL, x=NULL, title="Stage 9")+
#     scale_x_discrete(labels=c(expression("V"["e"]), expression("V"["g"]), expression("h"^2)))+
#     lims(y=c(0, .18))+
#     geom_errorbar(data=y[perm==0&(Source=="V(G)"|Source=="V(e)"|Source=="V(G)/Vp")], aes(x=Source, ymin=Variance-1.96*SE, ymax=Variance+1.96*SE), color="lightseagreen", width=0.05)
# 
# 
# 
# 
# 
# 
# z<-foreach(perm=c(0,1000:1999), .errorhandling="remove")%do%{
#     h<-fread(paste0("/mnt/pricey_2/priscilla/gcta/gcta_st7", perm, ".hsq"), fill=T)
#     h[,perm:=perm]
#     return(h)
# }
# 
# z<-rbindlist(z)
# 
# write.table(z,"/mnt/pricey_2/priscilla/gcta_hsq_1000perms_st7.txt", quote=F, sep="\t", row.names=F)

