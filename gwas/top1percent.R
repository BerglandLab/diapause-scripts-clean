
library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(foreach)
library(doMC)
registerDoMC(16)

#now read in GWAS, take top 1% with quantile rankings and return
files<-fread("/scratch/pae3g/genome-reconstruction/universal_input.txt")

y<-foreach(pop=files$V1[1:3000], phenotype=files$V2[1:3000],draw=files$V3[1:3000], perm=files$V4[1:3000], model=files$V6[1:3000], .errorhandling="remove")%dopar%{
    print(paste(draw, perm, sep=","))
    #read gwas adn do some fixes
    load(paste("/scratch/pae3g/revisions/genesis_", phenotype, "_draw", draw, "_perm", perm, "_pop", pop, "_" , model, "_allsnpgrm_wolbachia_dropmissing.Rdat", sep=""))
    gwas<-assoc.results
    gwas[,maf:=pmin(freq, 1-freq)]
    gwas<-gwas[maf>=0.05]
    gwas[,q:=log10(frank(Score.pval)/(length(Score.pval)+1))]
    gwas[,pop:=pop]
    gwas[,draw:=draw]
    return(gwas[q<=(-2)])
}

y<-rbindlist(y)
save(y, file="/scratch/pae3g/revisions/gwas_top1percent_dropmissing.Rdat")