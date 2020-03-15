#population genetics in lasso snps


library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(foreach)
library(doMC)
registerDoMC(10)
library(Rmisc)
library(ggbeeswarm)

files<-fread("/scratch/pae3g/genome-reconstruction/universal_input.txt")


p<-foreach(pop=files$V1[1:3000], phenotype=files$V2[1:3000],draw=files$V3[1:3000], perm=files$V4[1:3000], model=files$V6[1:3000], .errorhandling="remove")%dopar%{
    load(paste0("/scratch/pae3g/revisions/lasso/", "lasso_", phenotype, "_draw",draw, "_perm", perm, "_", model, "_pop", pop, ".lassoSites.Rdata"))
    return(sites)
}

p<-rbindlist(p)

p[, chr:=tstrsplit(site,split="_")[[1]]]
p[, pos:=as.numeric(tstrsplit(site,split="_")[[2]])]



#merge lasso.all with ZI data

zi<-fread("/scratch/pae3g/revisions/evolution/fst_ZI_AT_gr_12_fall.txt")

lasso.zi<-merge(p, zi, by=c("chr", "pos"))

lasso.zi[,pro.diapause.zambia:=ifelse(sign(coef)==1, p2, 1-p2)]

lasso.zi.sum<-lasso.zi[,.(min=min(pro.diapause.zambia, na.rm=T),
            med=median(pro.diapause.zambia, na.rm=T), 
            max=max(pro.diapause.zambia, na.rm=T),
            n=.N,
            prop.01=sum(pro.diapause.zambia>.01)/.N,
            prop.05=sum(pro.diapause.zambia>.05)/.N,
            prop.1=sum(pro.diapause.zambia>.1)/.N,
            prop.2=sum(pro.diapause.zambia>.2)/.N), .(pop, pheno, GRM, perm)]


lasso.zi.sum.melt<-melt(lasso.zi.sum, id.vars=c("pop", "pheno", "GRM", "perm"), measure.vars=c("med", "prop.01", "prop.05", "prop.1", "prop.2"))

lasso.zi.sum.melt[perm!=0&pop=="both", group:="both-permuted"]
lasso.zi.sum.melt[perm==0&pop=="both", group:="both-observed"]
lasso.zi.sum.melt[perm==0&pop=="A", group:="A-observed"]
lasso.zi.sum.melt[perm==0&pop=="B", group:="B-observed"]
lasso.zi.sum.melt[,pheno2:=ifelse(pheno=="diapause.bin", "stage 7", "stage 9")]

zi.plot<-ggplot(lasso.zi.sum.melt, aes(x=variable, y=value, color=group))+
    geom_beeswarm(cex=.25, dodge.width=.8,alpha=0.25)+
    labs(x="", y="freq. or proportion", color="")+
    scale_color_manual(values=c("#39568CFF", "#440154FF", "lightseagreen", "grey80"))+
    facet_grid(.~pheno2)+
    scale_x_discrete(labels=c("median", "> 0.01", "> 0.05", "> 0.10", "> 0.20"))+
    theme(strip.background = element_blank())


#test for enrichment of lasso snps in admixture tracts



load("/scratch/pae3g/revisions/evolution/admix_universal_lasso.Rdata")


admix.sum<-admix[,.(total.admix=sum(n.admix)/.N, unique.admix=sum(n.admix>0)/.N, n=.N), .(perm, draw, pheno, model, pop)]

admix.sum[perm!=0&pop=="both", group:="both-permuted"]
admix.sum[perm==0&pop=="both", group:="both-observed"]
admix.sum[perm==0&pop=="A", group:="A-observed"]
admix.sum[perm==0&pop=="B", group:="B-observed"]

admix.plot<-ggplot(admix.sum, aes(x=pheno, y=total.admix, color=group))+
    geom_beeswarm(cex=.25, dodge.width=.8,alpha=0.25)+
    labs(x="", y="admixture count per SNP", color="")+
    scale_color_manual(values=c("#39568CFF", "#440154FF", "lightseagreen", "grey80"))+
    scale_x_discrete(labels=c("stage 7", "stage 9"))+
    theme(legend.position="none")


a.7<-ecdf(admix.sum[group=="both-permuted"&pheno=="diapause.bin", unique.admix])
a.7(median(admix.sum[group=="both-observed"&pheno=="diapause.bin", unique.admix]))

a.9<-ecdf(admix.sum[group=="both-permuted"&pheno=="diapause.bin9", unique.admix])
a.9(median(admix.sum[group=="both-observed"&pheno=="diapause.bin9", unique.admix]))






#read in files


load("/scratch/pae3g/revisions/evolution/ihs_dgrp_universal_lasso.Rdata")

ihs.dgrp<-copy(p)
load("/scratch/pae3g/revisions/evolution/ihs_northern_universal_lasso.Rdata")

ihs.north<-copy(p)

a<-rbind(ihs.dgrp, ihs.north)
lasso.ihs.sum<-a[,.(min=min(IHS, na.rm=T),
                          med=median(IHS, na.rm=T), 
                          max=max(IHS, na.rm=T),
                          n=.N), .(pop, pheno, draw, perm, test, model)]


lasso.ihs.sum[perm!=0&pop=="both", group:="both-permuted"]
lasso.ihs.sum[perm==0&pop=="both", group:="both-observed"]
lasso.ihs.sum[perm==0&pop=="A", group:="A-observed"]
lasso.ihs.sum[perm==0&pop=="B", group:="B-observed"]

lasso.ihs.sum.melt<-melt(lasso.ihs.sum, measure.vars=c("min", "med", "max"), id.vars=c("draw", "pop", 'pheno', "test", "perm", "model", "group"))

lasso.ihs.sum.melt[,pheno2:=ifelse(pheno=="diapause.bin", "stage 7", "stage 9")]


ihs.plot<-ggplot(lasso.ihs.sum.melt, aes(x=variable, y=value, color=group))+
    geom_beeswarm(cex=.25, dodge.width=.8,alpha=0.25)+
    labs(x="", y="IHS", color="")+
    scale_color_manual(values=c("#39568CFF", "#440154FF", "lightseagreen", "grey80"))+
    facet_grid(test~pheno2)+
    scale_x_discrete(labels=c("minimum", "median", "maximum"))+ 
    theme(strip.text.x = element_blank(),strip.background = element_blank())



j<-foreach(testname=c("dgrp", "north", "dgrp", "north"), phenoname=c("diapause.bin", "diapause.bin", "diapause.bin9", "diapause.bin9" ))%do%{
    min.ecdf<-ecdf(lasso.ihs.sum[test==testname & pheno==phenoname & group=="both-permuted", min])
    min.q<-min.ecdf(median(lasso.ihs.sum[test==testname & pheno==phenoname&group=="both-observed", min]))
    med.ecdf<-ecdf(lasso.ihs.sum[test==testname&pheno==phenoname&group=="both-permuted", med])
    med.q<-med.ecdf(median(lasso.ihs.sum[test==testname&pheno==phenoname&group=="both-observed", med]))
    max.ecdf<-ecdf(lasso.ihs.sum[test==testname&pheno==phenoname&group=="both-permuted", max])
    max.q<-max.ecdf(median(lasso.ihs.sum[test==testname&pheno==phenoname&group=="both-observed", max]))
    return(data.table(test=testname,
                      pheno=phenoname,
                      min.q=min.q,
                      med.q=med.q,
                      max.q=max.q))
}
j<-(rbindlist(j))

