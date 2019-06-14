
# library(data.table)
# library(foreach)
# library(cowplot)
# library(doMC)
# registerDoMC(10)
# load(file="/mnt/spicy_1/pigmentation/inputData/dat.Rdata")
# setkey(dat, chr, pos, pop)
# ids=as.data.table(read.csv("/mnt/pricey_2/priscilla/hybrid/all_popinfo_formap_set.csv"))
# setnames(ids, "name", "pop")
# setkey(ids, pop)
# 
# dat<-merge(dat, ids, by="pop")
# 
# #subset to core20
# core=dat[set=="Core20"]
# 
# save(core, file="/mnt/pricey_2/priscilla/core20.rdat")
# #calculate ref freq
# core[,rf:=1-af]
# #logit transform ref freq
# core[,logit.rf:=qlogis(rf)]
# core[,pop.char:=as.character(pop_name)]
# 
# setkey(core, pop.char)
# deltas<-foreach(p=unique(core$pop.char))%do%{
#     print(p)
#     #calculate actual change and logit -transformed change in allele frequency. spring-fall so will be positive when spring freqeuncy is higher (same sign as gwas)
#     delta=core[pop.char==p,.(diff.logit=logit.rf[season=="spring"]-logit.rf[season=="fall"], diff=rf[season=="spring"]-rf[season=="fall"],population=p), .( chr, pos )]
#     return(delta)
# }
# names(deltas)=unique(core$pop.char)
# save(deltas, file= "/mnt/pricey_2/priscilla/core20delta.rdat")

#move to rivanna and code there now


library(data.table)
library(foreach)
library(cowplot)
library(doMC)
registerDoMC(20)

load("/scratch/pae3g/evolution/core20delta.rdat")

pops<-names(deltas)

pop.test<-foreach(pop=pops)%do%{
    
    clump.enrich<-foreach(perm =c(0,101:200)) %dopar% { #cycle through original data and 100 permutations
        print(perm)
        
        gwas<-fread(paste0("/scratch/pae3g/evolution/gwas_p_score_inv_id_perm", perm, ".txt"))
        seas<-deltas[[pop]]
        seas<-seas[diff.logit!=Inf&diff.logit!=(-Inf)]
        
        gwas<-merge(gwas, seas, by=c("chr", "pos"))
        
        #quantile rank seasonal changes
        gwas[,logit.q:=frank(diff.logit)/(length(diff.logit)+1)]
        gwas[,diff.q:=frank(diff)/(length(diff)+1)]
        
        #Z tranform seasonal changes
        gwas[,logit.Z.q:=qnorm(logit.q, 0,1)]
        gwas[,diff.Z.q:=qnorm(diff.q, 0,1)]
        
        #quantile rank and normalize Scores 
        gwas[,Score.q:=frank(gwas.Score)/(length(gwas.Score)+1)]
        gwas[,Score.Z.q:=qnorm(Score.q, 0,1)]
        
        #take product of Score * seasonal difference
        gwas[,Score.Z.diff.Z:=Score.Z.q*diff.Z.q]
        gwas[,Score.Z.diff:=Score.Z.q*diff]
        gwas[,Score.Z.logit.Z:=Score.Z.q*logit.Z.q]
        gwas[,Score.Z.logit:=Score.Z.q*diff.logit]
        
        #test for concordance of cline:
        gwas[,conc:=sign(diff)==sign(gwas.Score)]
        #take product of Score.Z*seas.beta
        a<-foreach (top = c(101, 490, 5316, 25250, 67914)) %do% { #these numbers correspond to the top # of snps in our different FDR cutoffs
            snps<-fread(paste0("/nv/vol186/bergland-lab/Priscilla/perm_summaries/perm", perm, "_top", top, "_ld0.4_200kb.clumped"))
            setnames(snps, c("CHR", "BP"), c("chr", "pos"))
            snps[,chr:=as.character(chr)]
            snps[chr=="23", chr:="X"]
            gwas.snps<-merge(gwas, snps[, .(chr, pos)], by=c("chr", "pos"))
            
            #return sums
            return(data.table(perm=perm,
                              top=top,
                              pop=pop,
                              sum.diff.Z=sum(gwas.snps$Score.Z.diff.Z),
                              sum.diff=sum(gwas.snps$Score.Z.diff),
                              sum.logit.Z.=sum(gwas.snps$Score.Z.logit.Z),
                              sum.logit=sum(gwas.snps$Score.Z.logit),
                              n=nrow(gwas.snps),
                              n.conc=sum(gwas.snps$conc, na.rm=T)))
            
        }
        return(rbindlist(a))
    }
    
    clump.enrich<-rbindlist(clump.enrich)
    clump.enrich[,prop.conc:=n.conc/n]
    return(clump.enrich)
}
 

pop.test<-rbindlist(pop.test)  
write.table(pop.test, "/scratch/pae3g/evolution/index_snps_PRS_bypop_seas2019.txt", quote=F, sep="\t", row.names=F)


library(data.table)
library(foreach)
library(cowplot)
library(ggbeeswarm)

c<-fread("/scratch/pae3g/evolution/index_snps_PRS_bypop_seas2019.txt")

c[,avg.diff:=sum.diff/n]
c[,avg.diff.Z:=sum.diff.Z/n]
c[,avg.logit:=sum.logit/n]
c[,avg.logit.Z:=sum.logit.Z./n]

c[,rank.avg.diff:=frank(avg.diff), .(top, pop)]
c[,rank.avg.diff.Z:=frank(avg.diff.Z), .(top, pop)]
c[,rank.avg.logit:=frank(avg.logit), .(top, pop)]
c[,rank.avg.logit.Z:=frank(avg.logit.Z), .(top, pop)]

c[perm==0&rank.avg.diff>96]
c[perm==0&rank.avg.diff.Z>96]
c[perm==0&rank.avg.logit>96]
c[perm==0&rank.avg.logit.Z>96]

c[perm==0&rank.avg.diff<10]
c[perm==0&rank.avg.diff.Z<10]
c[perm==0&rank.avg.logit<10]
c[perm==0&rank.avg.logit.Z<10]

ggplot(c[perm!=0], aes(x=as.factor(top), y=prop.conc))+
    geom_quasirandom(method="smiley", color="grey80", size=0.5)+
    geom_point(data=c[perm==0], aes(x=as.factor(top), y=prop.conc), color="lightseagreen", size=1)+
    scale_x_discrete(labels=c("FDR < 0.025", "FDR < 0.05", "FDR < 0.1", "FDR < 0.15", "FDR < 0.2"))+
    theme(axis.text.x=element_text(angle=45, hjust=1))+
    labs(x="", y="proportion concordant")+
    facet_wrap(~pop)

ggplot(c[perm!=0&top>500], aes(x=as.factor(top), y=avg.diff))+
    geom_quasirandom(method="smiley", color="grey80", size=0.5)+
    geom_point(data=c[perm==0&top>500], aes(x=as.factor(top), y=avg.diff), color="lightseagreen", size=1)+
    scale_x_discrete(labels=c("FDR < 0.025", "FDR < 0.05", "FDR < 0.1", "FDR < 0.15", "FDR < 0.2"))+
    theme(axis.text.x=element_text(angle=45, hjust=1))+
    labs(x="", y="avg seasonal score")+
    facet_wrap(~pop)

ggplot(c[perm!=0&top>500], aes(x=as.factor(top), y=avg.diff.Z))+
    geom_quasirandom(method="smiley", color="grey80", size=0.5)+
    geom_point(data=c[perm==0&top>500], aes(x=as.factor(top), y=avg.diff.Z), color="lightseagreen", size=1)+
    scale_x_discrete(labels=c("FDR < 0.025", "FDR < 0.05", "FDR < 0.1", "FDR < 0.15", "FDR < 0.2"))+
    theme(axis.text.x=element_text(angle=45, hjust=1))+
    labs(x="", y="avg seasonal score")+
    facet_wrap(~pop)


ggplot(c[perm!=0], aes(x=as.factor(top), y=avg.logit.Z))+
    geom_quasirandom(method="smiley", color="grey80", size=0.5)+
    geom_point(data=c[perm==0], aes(x=as.factor(top), y=avg.logit.Z), color="lightseagreen", size=1)+
    scale_x_discrete(labels=c("FDR < 0.025", "FDR < 0.05", "FDR < 0.1", "FDR < 0.15", "FDR < 0.2"))+
    theme(axis.text.x=element_text(angle=45, hjust=1))+
    labs(x="", y="avg seasonal score")+
    facet_wrap(~pop)+
    geom_hline(yintercept=0, color="lightseagreen", size=0.5)

ggplot(c[perm!=0&top>500], aes(x=as.factor(top), y=avg.logit))+
    geom_quasirandom(method="smiley", color="grey80", size=0.5)+
    geom_point(data=c[perm==0&top>500], aes(x=as.factor(top), y=avg.logit), color="lightseagreen", size=1)+
    scale_x_discrete(labels=c("FDR < 0.025", "FDR < 0.05", "FDR < 0.1", "FDR < 0.15", "FDR < 0.2"))+
    theme(axis.text.x=element_text(angle=45, hjust=1))+
    labs(x="", y="avg seasonal score")+
    facet_wrap(~pop)
