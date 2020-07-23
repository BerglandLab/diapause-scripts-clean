library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(foreach)
library(doMC)
registerDoMC(16)

#now read in GWAS, take top 1% with quantile rankings and return
files<-fread("/scratch/pae3g/genome-reconstruction/universal_input.txt")

y<-foreach(pop=files$V1[1:3000], phenotype=files$V2[1:3000],draw=files$V3[1:3000], perm=files$V4[1:3000], .errorhandling="remove")%dopar%{
    print(paste(draw, perm, sep=","))
    #read gwas adn do some fixes
    load(paste("/scratch/pae3g/revisions/genesis_", phenotype, "_draw", draw, "_perm", perm, "_pop", pop, "_loco_allsnpgrm_wolbachia_dropmissing.Rdat", sep=""))
    loco<-copy(assoc.results)
    load(paste("/scratch/pae3g/revisions/genesis_", phenotype, "_draw", draw, "_perm", perm, "_pop", pop, "_nonloco_allsnpgrm_wolbachia_dropmissing.Rdat", sep=""))
    nonloco<-copy(assoc.results)
    loco[,maf:=pmin(freq, 1-freq)]
    loco<-loco[maf>=0.05]
    nonloco[,maf:=pmin(freq, 1-freq)]
    nonloco<-nonloco[maf>=0.05]
    b<-merge(loco[,.(chr, pos, Score.pval)], nonloco[,.(chr, pos, Score.pval)], by=c("chr", "pos"), suffix=c(".loco", ".nonloco"))
    a<-cor(b$Score.pval.loco, b$Score.pval.nonloco, method="spearman")
    return(data.table(corr=a,
                      pop=pop,
                      perm=perm,
                      draw=draw,
                      pheno=phenotype))
}

y<-rbindlist(y)

save(y, file="/scratch/pae3g/revisions/loco-nonloco-corr.Rdat")


load("/scratch/pae3g/revisions/loco-nonloco-corr.Rdat")
y[,permuted:=ifelse(perm==0, F, T)]
         
corr.sum<-y[,.(n=.N, med.corr=median(corr), sd.corr=sd(corr)), .(permuted)]