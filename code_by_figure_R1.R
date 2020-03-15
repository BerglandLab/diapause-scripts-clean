#Erickson et al 2020 figures


#################################
#### FIGURE 1 ###################
#################################

#make a map of starting lines
library(maps)
library(data.table)
library(cowplot)

#these data are available in File S2

pops<-fread("~/Box Sync/HSparents/strain_geography.csv") #on dryad
geo.sum<-pops[,.(n=.N), .(latitude, longitude, geography)]
#move pa fall to the left a little bit so two points show up
geo.sum[geography=="Penn-fall", longitude:=-75.88]
geo.sum[,geography:=factor(geography, levels=c("Carribean", "Southeast", "North Carolina", "Penn-fall", "Penn-spring", "New York", "Maine"))]
world <- map_data("world")

ggplot()+ geom_polygon( data=world, aes(x=long, y=lat, group = group), colour="white", fill="lightgrey")+geom_point(data=geo.sum, aes(x=longitude, y=latitude, color=geography, size=as.numeric(n)) ,alpha=0.8)+labs(color="Origin of founder line", x="Longitude", y="Latitude", size="Number of lines")+coord_cartesian(xlim=c(-90,-65), ylim=c(20,45))+scale_size_continuous(breaks=c(1,2,10))

##############################
#### FIGURE 2 ################
##############################

library(data.table)
library(cowplot)
library(magick)
#on laptop: basic phenotypes

#read in phenotypes adn create categorical phenotypes based on ovary + eggs

p<-fread("~/Box Sync/hybridSwarm/phenos_wolbachia.txt") #on dryad
p[Ovary<=7&Eggs==0, diapause.group:="<=7"]
p[Ovary>7&Ovary<=9&Eggs==0, diapause.group:="8-9"]
p[Ovary==10&Eggs==0, diapause.group:="10"]
p[Ovary>10&Eggs==0, diapause.group:="11-14"]
p[Eggs==1, diapause.group:="1 Egg"]
p[Eggs>1&Eggs<=5, diapause.group:="2-5 Eggs"]
p[Eggs>5&Eggs<=10, diapause.group:="6-10 Eggs"]
p[Eggs>10&Eggs<15, diapause.group:="10-14 Eggs"]
p[Eggs>=15, diapause.group:="15+ Eggs"]

#get phenotypes for most advanced ovary stage in flies with eggs
p[Ovary<=7&Eggs>0, egg.group:="<=7"]

p[Ovary>7 & Ovary <= 9 & Eggs > 0, egg.group:="8-9"]
p[Ovary==10&Eggs>0, egg.group:="10"]
p[Ovary>10&Eggs>0, egg.group:="11-14"]

#reorder factors
p[,diapause.group:=factor(p$diapause.group, levels=c("<=7", "8-9", "10", "11-14", "1 Egg", "2-5 Eggs", "6-10 Eggs", "10-14 Eggs", "15+ Eggs"))]

p[,egg.group:=factor(p$egg.group, levels=c("<=7", "8-9", "10", "11-14", "1 Egg", "2-5 Eggs", "6-10 Eggs", "10-14 Eggs", "15+ Eggs"))]



#group temperatures
p[,temp.group:=NULL]
p[temp.rack.cal<11, temp.group:="10-11°C"]
p[temp.rack.cal>=11 & temp.rack.cal<12, temp.group:="11-12°C"]
p[temp.rack.cal>=12 & temp.rack.cal<13, temp.group:="12-13°C"]
p[temp.rack.cal>=13 & temp.rack.cal<14, temp.group:="13-14°C"]
p[temp.rack.cal>=14, temp.group:="14°C +"]


a<-ggplot()
#plot ovary phenotype as a function of temperature group
b.sum<-p[!is.na(diapause.group), .(n=.N), .(temp.group)]
b.sum[,y:=1.02]
b<-ggplot(p[!is.na(diapause.group)], aes(x=temp.group, fill=diapause.group))+geom_bar(position="fill")+
         labs(x="Temperature", 
         fill="Ovary\nPhenotype", 
         y="Proportion of Flies", 
         title="All flies")+
    theme(axis.text.x=element_text(angle=45,hjust=1))
#plot ovary stages of individuals with eggs
c.sum<-p[!is.na(egg.group), .(n=.N), .(temp.group)]
c<-ggplot(p[!is.na(egg.group)], aes(x=temp.group, fill=egg.group))+geom_bar(position="fill")+
    labs(x="Temperature", 
         fill="Ovary\nPhenotype", 
         y="Proportion of Flies", 
         title="Flies with eggs")+
    theme(axis.text.x=element_text(angle=45,hjust=1))+
    scale_fill_discrete(drop=F,limits = levels(p$egg.group), breaks=c("<=7", "8-9", "10", "11-14"))


#function for smoothed binomial plots
binomial_smooth <- function(...) {
    geom_smooth(method = "glm", method.args = list(family = "binomial"), se=FALSE)
}

#melt data across stages

p.melt<-melt(p, id.vars=c("id", "temp.rack.cal", "photoperiod"), measure.vars=c("diapause.bin9", "diapause.bin", "eggP"), variable.name="pheno", value.name="diapause")
p.melt[,pheno:=factor(p.melt$pheno, levels=c("eggP", "diapause.bin9", "diapause.bin"))]

#add jitter by photoperiod
p[photoperiod==9, diapause.plot:=as.numeric(diapause.bin)+.02]
p[photoperiod==11, diapause.plot:=as.numeric(diapause.bin)+.01]
p[photoperiod==13, diapause.plot:=as.numeric(diapause.bin)-.01]
p[photoperiod==15, diapause.plot:=as.numeric(diapause.bin)-.02]

p[photoperiod==9, diapause.bin9.plot:=as.numeric(diapause.bin9)+.02]
p[photoperiod==11, diapause.bin9.plot:=as.numeric(diapause.bin9)+.01]
p[photoperiod==13, diapause.bin9.plot:=as.numeric(diapause.bin9)-.01]
p[photoperiod==15, diapause.bin9.plot:=as.numeric(diapause.bin9)-.02]

p[photoperiod==9, egg.plot:=as.numeric(eggP)+.02]
p[photoperiod==11, egg.plot:=as.numeric(eggP)+.01]
p[photoperiod==13, egg.plot:=as.numeric(eggP)-.01]
p[photoperiod==15, egg.plot:=as.numeric(eggP)-.02]

#plot models of likelihood of diapause as a function of temperature for different photoperiods
d<-ggplot(p, aes(x=temp.rack.cal, y=diapause.bin, color=as.factor(photoperiod)))+
    binomial_smooth()+
    labs(x="Temperature °C", 
         y="Probability of diapause", 
         color="Photoperiod",
         title="Stage 8")+
    lims(y=c(-.05,1.05))+
    geom_point(data=p, aes(x=temp.rack.cal, y=diapause.plot, color=as.factor(photoperiod)), alpha=0.2)


e<-ggplot(p, aes(x=temp.rack.cal, y=diapause.bin9, color=as.factor(photoperiod)))+
    binomial_smooth()+
    labs(x="Temperature °C", 
         y="Probability of diapause", 
         color="Photoperiod",
         title="Stage 10")+
    lims(y=c(-.05,1.05))+
    geom_point(data=p, aes(x=temp.rack.cal, y=diapause.bin9.plot, color=as.factor(photoperiod)), alpha=0.2)

f<-ggplot(p, aes(x=temp.rack.cal, y=eggP, color=as.factor(photoperiod)))+
    binomial_smooth()+
    labs(x="Temperature °C",
         y="Probability of diapause", 
         color="Photoperiod",
         title="Stage 14")+
    lims(y=c(-.05,1.05))+
    geom_point(data=p, aes(x=temp.rack.cal, y=egg.plot, color=as.factor(photoperiod)), alpha=0.2)

leg<-get_legend(b)

#plot top row of multipanel plot

top.row<-plot_grid(a,
                   b+ theme(legend.position="none"),
                   c+ theme(legend.position="none"), 
                   leg,
                   nrow=1, labels=c("A", "B", "C"), rel_widths=c(0.6, 1.2, 1.2, .6))

#make bottom row
prow <- plot_grid( d + theme(legend.position="none"),
                   e + theme(legend.position="none"),
                   f + theme(legend.position="none"),
                   labels = c("D", "E", "F"),
                   nrow = 1
)


#add legend
legend <- get_legend(d)

bottom.row <- plot_grid( prow, legend, rel_widths = c(3, .6))

#FINAL PLOT:
pdf("~/Box Sync/manuscripts/diapause gwas/R1/figure2.pdf", height=6, width=8)
plot_grid(top.row, bottom.row, nrow=2)
dev.off()


#stats for figure legend for this and supplemental figures

summary(aov(glm(diapause.bin~temp.rack.cal+(photoperiod)+swarm+generation, data=p, family="binomial")))
summary(aov(glm(diapause.bin9~temp.rack.cal+(photoperiod)+swarm+generation, data=p, family="binomial")))
summary(aov(glm(eggP~temp.rack.cal+(photoperiod)+swarm+generation, data=p, family="binomial")))



#what is correlation between each phenotype?
p<-p[!is.na(diapause.bin)&!is.na(diapause.bin9)&!is.na(eggP)]
cor.test(p$diapause.bin, p$diapause.bin9)[["p.value"]]
cor.test(p$diapause.bin, p$eggP)[["p.value"]]
cor.test(p$diapause.bin9, p$eggP)[["p.value"]]

############################
######## FIGURE 3 #########
############################


library(foreach)
library(data.table)
library(ROCR)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(doMC)
registerDoMC(10)
library(ggbeeswarm)



y<-fread("/scratch/pae3g/revisions/gcta_hsq_1000perms.txt") #on dryad

h.plot<-ggplot(y[perm!=0&(Source=="V(G)"|Source=="V(e)"|Source=="V(G)/Vp")], aes(x=Source, y=Variance))+
    geom_quasirandom(color="grey80", varwidth = TRUE, size=.5, method = "smiley")+
    geom_point(data=y[perm==0&(Source=="V(G)"|Source=="V(e)"|Source=="V(G)/Vp")], aes(x=Source, y=Variance), color="lightseagreen", size=2)+
    labs(y="var. or prop.", x=NULL)+
    scale_x_discrete(labels=c(expression("V"["e"]), expression("V"["g"]), expression("h"^2)))+
    lims(y=c(0, .18))+
    geom_errorbar(data=y[perm==0&(Source=="V(G)"|Source=="V(e)"|Source=="V(G)/Vp")], aes(x=Source, ymin=Variance-1.96*SE, ymax=Variance+1.96*SE), color="lightseagreen", width=0.05)



z<-fread("/scratch/pae3g/revisions/gcta_hsq_1000perms_st7.txt") #on dryad

h7.plot<-ggplot(z[perm!=0&(Source=="V(G)"|Source=="V(e)"|Source=="V(G)/Vp")], aes(x=Source, y=Variance))+
    geom_quasirandom(color="grey80", varwidth = TRUE, size=.5, method = "smiley")+
    geom_point(data=z[perm==0&(Source=="V(G)"|Source=="V(e)"|Source=="V(G)/Vp")], aes(x=Source, y=Variance), color="lightseagreen", size=2)+
    labs(y="var. or prop.", x=NULL)+scale_x_discrete(labels=c(expression("V"["e"]), expression("V"["g"]), expression("h"^2)))+lims(y=c(0, .18))+
    geom_errorbar(data=z[perm==0&(Source=="V(G)"|Source=="V(e)"|Source=="V(G)/Vp")], aes(x=Source, ymin=Variance-1.96*SE, ymax=Variance+1.96*SE), color="lightseagreen", width=0.05)


#GIFs
files<-fread("/scratch/pae3g/genome-reconstruction/universal_input.txt") #on dryad

p<-foreach(pop=files$V1[1:3000], phenotype=files$V2[1:3000],draw=files$V3[1:3000], perm=files$V4[1:3000], model=files$V6[1:3000], .errorhandling="remove")%dopar%{
    load(paste0("/scratch/pae3g/revisions/lasso/", "GIF_", phenotype, "_draw",draw, "_perm", perm, "_", model, "_pop", pop, ".Rdata"))
    return(g)
}

p<-rbindlist(p)

save(p, file="/scratch/pae3g/revisions/GIF.Rdata") #on dryad

p[perm!=0&pop=="both", group:="both-permuted"]
p[perm==0&pop=="both", group:="both-observed"]
p[perm==0&pop=="A", group:="A-observed"]
p[perm==0&pop=="B", group:="B-observed"]

g.plot<-ggplot(p[!is.na(lambda.05)&!is.na(group)], aes(x=pheno, y=lambda.05, color=group))+
    geom_beeswarm(cex=.5, dodge.width = .5, alpha=0.5, size=0.5)+
    labs(y="λ GC", x=NULL, color=NULL)+
    scale_color_manual(values=c("dodgerblue2", "darkorchid", "lightseagreen", "grey80"))+
    #theme(legend.position = "none")+
    scale_x_discrete(labels=c("stage 8", "stage 10"))+
    lims(y=c(0.75, 1.15))
s


r1<-plot_grid(h7.plot, h.plot, g.plot, nrow=1, align="vh", axis="b", labels=c("A", "B", "C"), rel_widths=c(0.25, 0.25, 0.5))



o<-foreach(pop=files$V1[1:3000], phenotype=files$V2[1:3000],draw=files$V3[1:3000], perm=files$V4[1:3000], model=files$V6[1:3000], .errorhandling="remove")%dopar%{
    load(paste0("/scratch/pae3g/revisions/lasso/", "lasso_", phenotype, "_draw",draw, "_perm", perm, "_", model, "_pop", pop, ".lassoPred.Rdata"))
    return(pred)
}

o<-rbindlist(o)

save(o, file="/scratch/pae3g/revisions/lassoPred.Rdata") #on dryad

setkey(o, lassoMod)
o <- o[J(c("env", "env.pc", "env.pc.geno"))]



get_RORC_stat  <-  function(pred.i, obs.i, which.param="auc") {
    # pred.i <- pred[lassoMod=="env"]$pred; obs.i <- pred[lassoMod=="env"]$obs; which.param="auc"
    
    pred.obj <- prediction(pred.i, obs.i)
    perf.obj <- performance(pred.obj, which.param)
    
    return(perf.obj@y.values[[1]])
    
}


o.auc <- o[,list(auc=get_RORC_stat(pred.i=pred, obs.i=obs)), list(lassoMod, maf, perm,  loco, GRM, pheno, pop)]

o.auc[perm!=0&pop=="both", group:="both-permuted"]
o.auc[perm==0&pop=="both", group:="both-observed"]
o.auc[perm==0&pop=="A", group:="A-observed"]
o.auc[perm==0&pop=="B", group:="B-observed"]

o.roc <- o[,list(tpr=get_RORC_stat(pred.i=pred, obs.i=obs, "tpr"),
                 fpr=get_RORC_stat(pred.i=pred, obs.i=obs, "fpr")),
           list(lassoMod, maf, perm,  loco, GRM, pheno, pop)]


o.roc[perm!=0&pop=="both", group:="both-permuted"]
o.roc[perm==0&pop=="both", group:="both-observed"]
o.roc[perm==0&pop=="A", group:="A-observed"]
o.roc[perm==0&pop=="B", group:="B-observed"]

o.auc.plot <- ggplot(data=o.auc[lassoMod!="env"&pheno=="diapause.bin9"][order(perm, decreasing=T)], aes(x=lassoMod, y=auc, color=group)) +
    geom_beeswarm(cex=.25, dodge.width = .5, alpha=0.5, size=0.5) + 
    theme(legend.position = "none")+
    geom_beeswarm(data=o.auc[pheno=="diapause.bin9"&lassoMod=="env"&(perm==1001|(perm==0&GRM==1))], aes(x=lassoMod, y=auc, color=group),  cex=.25, dodge.width = .1, alpha=0.5, size=0.5)+
    scale_color_manual(values=c("dodgerblue2", "darkorchid", "lightseagreen", "grey80"))+
    labs(x=NULL, y="AUC")+scale_x_discrete(labels=c("env", "env + PC", "env + PC\n+ GWAS"))

o.roc.plot <- ggplot(data=o.roc[pheno=="diapause.bin9"&lassoMod=="env.pc.geno"][order(perm, decreasing=T)], aes(x=fpr, y=tpr, color=group, group=interaction(perm, GRM))) +
    geom_line(size=0.1) +
    theme(legend.position = "none")+
    labs(x="FPR", y="TPR")+
    scale_color_manual(values=c("dodgerblue2", "darkorchid", "lightseagreen", "grey80"))


p.auc.plot <- ggplot(data=o.auc[lassoMod!="env"&pheno=="diapause.bin"][order(perm, decreasing=T)], aes(x=lassoMod, y=auc, color=group)) +
    geom_beeswarm(cex=.25, dodge.width = .5, alpha=0.5, size=0.5) + 
    theme(legend.position = "none")+
    geom_beeswarm(data=o.auc[pheno=="diapause.bin"&lassoMod=="env"&(perm==1001|(perm==0&GRM==1))], aes(x=lassoMod, y=auc, color=group),  cex=.25, dodge.width = .1, alpha=0.5, size=0.5)+
    scale_color_manual(values=c("dodgerblue2", "darkorchid", "lightseagreen", "grey80"))+
    labs(x=NULL, y="AUC")+scale_x_discrete(labels=c("env", "env + PC", "env + PC\n+ GWAS"))

p.roc.plot <- ggplot(data=o.roc[pheno=="diapause.bin"&lassoMod=="env.pc.geno"][order(perm, decreasing=T)], aes(x=fpr, y=tpr, color=group, group=interaction(perm, GRM))) +
    geom_line(size=0.1) +
    theme(legend.position = "none")+
    labs(x="FPR", y="TPR")+
    scale_color_manual(values=c("dodgerblue2", "darkorchid", "lightseagreen", "grey80"))



#choose a random imputation to sample for manhattan plot

draw<-95
#(draw 95 chosen and looks very nice)

perm=0

load(paste("/scratch/pae3g/revisions/genesis_diapause.bin9_draw", draw, "_perm0_popboth_nonloco_allsnpgrm_wolbachia.Rdat", sep=""))
gwas<-assoc.results
gwas[,maf:=pmin(freq, 1-freq)]
gwas<-gwas[maf>=0.05]

gwas[,color:=1]
gwas[chr=="2R"|chr=="3R", color:=2]

#load in lasso snps for this gwas
load(paste0("/scratch/pae3g/revisions/lasso/lasso_diapause.bin9_draw",draw, "_perm0_nonloco_popboth.lassoSites.Rdata"))

sites[, chr:=tstrsplit(site,split="_")[[1]]]
sites[, pos:=as.numeric(tstrsplit(site,split="_")[[2]])]

sites<-merge(sites, gwas, by=c("chr", "pos"))

chr.table<-gwas[,.(max.id=max(variant.id)), .(chr)]
chr.table<-foreach(i=1:5, .combine="rbind")%do%{
    temp<-chr.table[i]
    temp[,midpoint:=(chr.table[i,max.id]+chr.table[i-1, max.id])/2]
    return(temp)
}
chr.table[chr=="2L", midpoint:=max.id/2]


mh<-ggplot(gwas, aes(x=variant.id, y=-log10(Score.pval), color=as.factor(color)))+
    geom_point()+
    scale_color_manual(values=c("black", "grey20"))+
    geom_point(data=sites,  aes(x=variant.id, y=-log10(Score.pval)), color="lightseagreen")+
    theme(legend.position="none")+labs(x=NULL, y=expression("-log"[10]*"(P)"))+scale_x_continuous(breaks=chr.table$midpoint, labels=c("2L", "2R", "3L", "3R", "X"))

r2<-mh

r3<-plot_grid(p.roc.plot, o.roc.plot, labels=c("E", "F"))

r4<-plot_grid(p.auc.plot, o.auc.plot, labels=c("G", "H"))


jpeg("/scratch/pae3g/revisions/figures/test_figure3.jpeg", height=8, width=8, res=300, units="in")
plot_grid(r1, r2, r3, r4, nrow=4, labels=c("", "D", "", ""))
dev.off()




####################################
##### FIGURE 4 #####################
####################################

#hs = hybrid swarm; d= dgrp  hp= hybrid parents
#r= rare allels maf <.05
#c = common alleles maf > .1
library(data.table)
library(Rmisc)

#hybrid swarm data
hs.c<-fread("/scratch/pae3g/evolution/ld_decay_mafover0.1_byswarm_revised.txt") #on dryad
hs.r<-fread("/scratch/pae3g/evolution/ld_decay_mafless0.05_byswarm_revised.txt") #on dryad

hs.c[,r2:=ld*ld]
hs.c[,r2a:=ld.a*ld.a]
hs.c[,r2b:=ld.b*ld.b]
hs.c.sum<-hs.c[!is.na(r2), .(med.r2=median(r2, na.rm=T),
                             mu.r2=mean(r2, na.rm=T),
                             ci.r2.upper=CI(r2, ci=.95)[1],
                             ci.r2.lower=CI(r2, ci=.95)[3],
                             med.r2a=median(r2a, na.rm=T),
                             med.r2b=median(r2b, na.rm=T)), .(dist)]

hs.r[,r2:=ld*ld]
hs.r[,r2a:=ld.a*ld.a]
hs.r[,r2b:=ld.b*ld.b]
hs.r.sum<-hs.r[!is.na(r2), .(med.r2=median(r2, na.rm=T),
                             mu.r2=mean(r2, na.rm=T),
                             ci.r2.upper=CI(r2, ci=.95)[1],
                             ci.r2.lower=CI(r2, ci=.95)[3],
                             med.r2a=median(r2a, na.rm=T),
                             med.r2b=median(r2b, na.rm=T)), .(dist)]

#dgrp data
d.r<-fread("/scratch/pae3g/evolution/ld_decay_dgrp_mafless0.05_revised.txt") #on dryad
d.r[,r2:=ld*ld]

d.r.sum<-d.r[!is.na(r2), .(med.r2=median(r2, na.rm=T), 
                           mu.r2=mean(r2, na.rm=T),
                           ci.r2.upper=CI(r2, ci=.95)[1],
                           ci.r2.lower=CI(r2, ci=.95)[3]), .(dist)]

d.c<-fread("/scratch/pae3g/evolution/ld_decay_dgrp_mafover0.1_revised.txt") #on dryad
d.c[,r2:=ld*ld]

d.c.sum<-d.c[!is.na(r2), .(med.r2=median(r2, na.rm=T),
                           mu.r2=mean(r2, na.rm=T),
                           ci.r2.upper=CI(r2, ci=.95)[1],
                           ci.r2.lower=CI(r2, ci=.95)[3]), .(dist)]

#downsampled dgrp data
d.r.d<-fread("/scratch/pae3g/evolution/ld_decay_dgrp_mafless0.05_downsampled_revised.txt") #on dryad
d.r.d[,r2:=ld*ld]

d.r.d.sum<-d.r.d[!is.na(r2), .(med.r2=median(r2, na.rm=T), 
                               mu.r2=mean(r2, na.rm=T),
                               ci.r2.upper=CI(r2, ci=.95)[1],
                               ci.r2.lower=CI(r2, ci=.95)[3]), .(dist)]

d.c.d<-fread("/scratch/pae3g/evolution/ld_decay_dgrp_mafover0.1_downsampled_revised.txt") #on dryad
d.c.d[,r2:=ld*ld]

d.c.d.sum<-d.c.d[!is.na(r2), .(med.r2=median(r2, na.rm=T),
                               mu.r2=mean(r2, na.rm=T),
                               ci.r2.upper=CI(r2, ci=.95)[1],
                               ci.r2.lower=CI(r2, ci=.95)[3]), .(dist)]

#hybrid swarm parents
hp.r<-fread("/scratch/pae3g/evolution/ld_decay_mafless0.05_parents_revised.txt") #on dryad
hp.r[,r2:=ld*ld]

hp.r.sum<-hp.r[!is.na(r2), .(med.r2=median(r2, na.rm=T),
                             mu.r2=mean(r2, na.rm=T),
                             ci.r2.upper=CI(r2, ci=.95)[1],
                             ci.r2.lower=CI(r2, ci=.95)[3]), .(dist)]

hp.c<-fread("/scratch/pae3g/evolution/ld_decay_mafover0.1_parents_revised.txt") #on dryad
hp.c[,r2:=ld*ld]

hp.c.sum<-hp.c[!is.na(r2), .(med.r2=median(r2, na.rm=T),
                             mu.r2=mean(r2, na.rm=T),
                             ci.r2.upper=CI(r2, ci=.95)[1],
                             ci.r2.lower=CI(r2, ci=.95)[3]), .(dist)]


a<-ggplot(hp.c.sum, aes(x=log10(dist), y=log10(mu.r2)))+
    geom_line(color="deeppink1")+
    geom_line(data=hs.c.sum, aes(x=log10(dist), y=log10(mu.r2)), color="lightseagreen" )+
    geom_line(data=d.c.sum, aes(x=log10(dist), y=log10(mu.r2)), color="darkorange2" )+
    geom_line(data=d.c.d.sum, aes(x=log10(dist), y=log10(mu.r2)), color="goldenrod2" )+
    
    geom_ribbon(data=hs.c.sum, aes(x=log10(dist), ymin=log10(ci.r2.lower), ymax=log10(ci.r2.upper)), fill="lightseagreen", alpha=0.2)+
    geom_ribbon(data=d.c.sum, aes(x=log10(dist), ymin=log10(ci.r2.lower), ymax=log10(ci.r2.upper)), fill="darkorange2", alpha=0.2)+
    geom_ribbon(data=hp.c.sum, aes(x=log10(dist), ymin=log10(ci.r2.lower), ymax=log10(ci.r2.upper)), fill="deeppink1", alpha=0.2)+
    geom_ribbon(data=d.c.d.sum, aes(x=log10(dist), ymin=log10(ci.r2.lower), ymax=log10(ci.r2.upper)), fill="goldenrod2", alpha=0.2)+
    labs(x=expression("log"[10]*"(distance (bp))"), y=expression("log"[10]*"(R"^2*")"))+
    lims(y=c(-2.5, 0))


b<-ggplot(hp.r.sum, aes(x=log10(dist), y=log10(mu.r2)))+
    geom_line(color="deeppink1")+
    geom_line(data=hs.r.sum, aes(x=log10(dist), y=log10(mu.r2)), color="lightseagreen" )+
    geom_line(data=d.r.sum, aes(x=log10(dist), y=log10(mu.r2)), color="darkorange2" )+
    geom_line(data=d.r.d.sum, aes(x=log10(dist), y=log10(mu.r2)), color="goldenrod2" )+
    
    geom_ribbon(data=hs.r.sum, aes(x=log10(dist), ymin=log10(ci.r2.lower), ymax=log10(ci.r2.upper)), fill="lightseagreen", alpha=0.2)+
    geom_ribbon(data=d.r.sum, aes(x=log10(dist), ymin=log10(ci.r2.lower), ymax=log10(ci.r2.upper)), fill="darkorange2", alpha=0.2)+
    geom_ribbon(data=hp.r.sum, aes(x=log10(dist), ymin=log10(ci.r2.lower), ymax=log10(ci.r2.upper)), fill="deeppink1", alpha=0.2)+
    geom_ribbon(data=d.r.d.sum, aes(x=log10(dist), ymin=log10(ci.r2.lower), ymax=log10(ci.r2.upper)), fill="goldenrod2", alpha=0.2)+
    labs(x=expression("log"[10]*"(distance (bp))"), y=expression("log"[10]*"(R"^2*")"))+
    lims(y=c(-2.5, 0))

leftside<-plot_grid(a,b, labels=c("A", "B"), nrow=2, align="v", axis="l")
#hybrid swarm long distance LD

inter.hs=fread("/scratch/pae3g/revisions/interchromosomal_ld.txt") #on dryad
intra.hs=fread("/scratch/pae3g/revisions/chromsome_arm_ld.txt") #on dryad

intra.p=fread("/scratch/pae3g/revisions/chromsome_arm_ld_parents.txt") #on dryad
inter.p=fread("/scratch/pae3g/revisions/interchromosomal_ld_parents.txt") #on dryad

intra.d=fread("/scratch/pae3g/revisions/chromsome_arm_ld_dgrp.txt") #on dryad
inter.d=fread("/scratch/pae3g/revisions/interchromosomal_ld_dgrp.txt") #on dryad

inter.hs[,pop:="hs"]
intra.hs[,pop:="hs"]
inter.p[,pop:="p"]
intra.p[,pop:="p"]
inter.d[,pop:="d"]
intra.d[,pop:="d"]

intra<-rbind(intra.hs, intra.p, intra.d)
inter<-rbind(inter.hs, inter.p, inter.d)


c<-ggplot(data=intra[chr==2], aes(x=pop, y=log10(ld*ld), color=pop ))+
    geom_boxplot()+
    labs(x="", y=expression("log"[10]*"(R"^2*")"))+
    scale_color_manual(values=c("darkorange2", "lightseagreen", "deeppink1"))+
    scale_x_discrete(labels=c("DGRP", "hybrids", "parents"))+ 
    theme(legend.position = "none")+
    lims(y=c(-11,0))


d<-ggplot(data=intra[chr==3], aes(x=pop, y=log10(ld*ld), color=pop ))+
    geom_boxplot()+
    labs(x="", y=expression("log"[10]*"(R"^2*")"))+
    scale_color_manual(values=c("darkorange2", "lightseagreen", "deeppink1"))+
    scale_x_discrete(labels=c("DGRP", "hybrids", "parents"))+ 
    theme(legend.position = "none")+
    lims(y=c(-11,0))


e<-ggplot(data=inter, aes(x=pop, y=log10(ld*ld), color=pop ))+
    geom_boxplot()+
    labs(x="", y=expression("log"[10]*"(R"^2*")"))+
    scale_color_manual(values=c("darkorange2", "lightseagreen", "deeppink1"))+
    scale_x_discrete(labels=c("DGRP", "hybrids", "parents"))+ 
    theme(legend.position = "none")+
    lims(y=c(-11,0))


rightside<-plot_grid(c,d,e, nrow=3, labels=c("C", "D", "E"), align="v", axis="l")


pdf("/scratch/pae3g/revisions/figures/ld_decay.pdf", heigh=6, width=6)
plot_grid(leftside,rightside, nrow=1)
dev.off()


###################
### FIGURE 5 ######
###################


load("/scratch/pae3g/revisions/evolution/bergland2014_sign_universal_threshold.Rdata") #on dryad
#loads in as y
load("/scratch/pae3g/oldscratch_recovered/evolution/6d_data.Rdata")
b<-as.data.table(p)

load("/scratch/pae3g/revisions/all_lasso_for_alan.Rdata") #on dryad
#loads in as p

b[,maf:=pmin(f.hat, 1-f.hat)]
b<-b[maf>=0.05]
b<-b[clinal.beta<3&clinal.beta>(-3)]
b<-b[sfsfsfX.beta<3&sfsfsfX.beta>(-3)]

x<-merge(p, b[,.(chr, pos, clinal.beta,  maf, sfsfsfX.p, sfsfsfX.beta)], by=c("chr", "pos"))

x[,ps.cline:=-1*coef*clinal.beta]
x[,ps.seas:=1*coef*sfsfsfX.beta]

x.sum<-x[,.(clinal=sum(ps.cline),
            seasonal=sum(ps.seas)),. (GRM, perm, pop, pheno)]

x.melt<-melt(x.sum, measure.vars=c("clinal", "seasonal"), id.vars=c("GRM", "perm", "pop", "pheno"))
setnames(x.melt, c("variable", "value", "GRM"), c("test", "poly" , "draw" ))
x.melt[,top:="lasso"]

b2014<-rbind(y[,.(poly, top, perm, draw, test, pheno, pop)], x.melt)
b2014[perm!=0&pop=="both", group:="both-permuted"]
b2014[perm==0&pop=="both", group:="both-observed"]
b2014[perm==0&pop=="A", group:="A-observed"]
b2014[perm==0&pop=="B", group:="B-observed"]
b2014[perm!=0&pop=="A", group:="A-permuted"]
b2014[perm!=0&pop=="B", group:="B-permuted"]

b2014[,pheno2:=ifelse(pheno=="diapause.bin", "st. 8", "st. 10")]
b2014[,pheno2:=factor(b2014$pheno2, levels=c("st. 8", "st. 10"))]

b2014.sum<-b2014[,.(med=median(poly), q.025=quantile(poly, 0.025), q.975=quantile(poly, .975)), .(group, pheno, pheno2, pop, test,top)]
b2014.sum[,pheno2:=factor(b2014.sum$pheno2, levels=c("st. 8", "st. 10"))]

b2014.sum[top==-5, th:="Top 0.001%"]
b2014.sum[top==-4, th:="Top 0.01%"]
b2014.sum[top==-3, th:="Top 0.1%"]
b2014.sum[top==-2, th:="Top 1%"]
b2014.sum[top==-1, th:="Top 10%"]
b2014.sum[top==0, th:="all SNPs"]
b2014.sum[top=="lasso", th:="LASSO"]

a.plot<-ggplot(b2014.sum[test=="clinal" & th%in%c("LASSO", "all SNPs", "Top 0.01%", "Top 0.1%")])+
    geom_point(data=b2014.sum[test=="clinal"&th%in%c("LASSO", "all SNPs", "Top 0.01%", "Top 0.1%")], aes(x=pop, y=med, color=group), position=position_dodge(width=0.5))+
    geom_errorbar(data=b2014.sum[test=="clinal"&th%in%c("LASSO", "all SNPs", "Top 0.01%", "Top 0.1%")], aes(x=pop, ymax=q.975, ymin=q.025, color=group), width=0.2, position=position_dodge(width=0.5))+
    labs(x="", y="sum(GWAS coefficient*model coefficient)", color="", title="clinal")+
    scale_color_manual(values=c("dodgerblue2", "grey80", "darkorchid", "grey80", "lightseagreen", "grey80"))+
    facet_grid(th~pheno2, scales ="free_y")+
    theme(legend.position = "none")

b.plot<-ggplot(b2014.sum[test=="seasonal" & th%in%c("LASSO", "all SNPs", "Top 0.01%", "Top 0.1%")])+
    geom_point(data=b2014.sum[test=="seasonal"&th%in%c("LASSO", "all SNPs", "Top 0.01%", "Top 0.1%")], aes(x=pop, y=med, color=group), position=position_dodge(width=0.5))+
    geom_errorbar(data=b2014.sum[test=="seasonal"&th%in%c("LASSO", "all SNPs", "Top 0.01%", "Top 0.1%")], aes(x=pop, ymax=q.975, ymin=q.025, color=group), width=0.2, position=position_dodge(width=0.5))+
    labs(x="", y="", color="", title="seasonal")+
    scale_color_manual(values=c("dodgerblue2", "grey80", "darkorchid", "grey80", "lightseagreen", "grey80"))+
    facet_grid(th~pheno2, scales ="free_y")+
    theme(legend.position = "none")

pdf("/scratch/pae3g/revisions/figures/bergland2014_main.pdf", height=8, width=8)
plot_grid(a.plot, b.plot, labels=c("A", "B"), nrow=1)
dev.off()



###############################
#### FIGURE 6 #################
###############################



library(data.table)
library(cowplot)
library(lme4)

f0<-fread("~/Box Sync//manuscripts/diapause gwas/dryad_upload/mesocosm_F0.csv") #on dryad
f0[,collection.date:=mdy(Date)]
f0[Eggs==0&Ovary<10, diapause.bin9:=1]
f0[Eggs>0|(Eggs==0&Ovary>=10), diapause.bin9:=0]
f0[Eggs==0&Ovary<8, diapause.bin:=1]
f0[Eggs>0|(Eggs==0&Ovary>=8), diapause.bin:=0]
f0<-f0[!is.na(collection.date)&!is.na(Cage)]
f0[Cage%in%c(1:6), food:="fruit"]
f0[Cage%in%c(7:11), food:="cornmeal-mollasses"]

f0.sum<-f0[,.(prop.d9=sum(diapause.bin9, na.rm=T)/.N, 
              prop.d7=sum(diapause.bin, na.rm=T)/.N, 
              med.eggs=as.double(median(Eggs, na.rm=T)), 
              q75Eggs=quantile(Eggs, .75, na.rm=T),
              mad.eggs=mad(x=Eggs, constant=1/quantile(Eggs,.75, na.rm=T)),
              n=.N), .(collection.date, food)]
#add SE
f0.sum[,se.prop.d9:=sqrt(prop.d9*(1-prop.d9)/n)]
f0.sum[,se.prop.d7:=sqrt(prop.d7*(1-prop.d7)/n)]



a<-ggplot(f0.sum[collection.date!="2019-01-25"&food=="fruit"], aes(x=collection.date, y=prop.d9, color="fruit-reared\noutdoor F0s\nstage 10")) +
    geom_point()+
    geom_line()+
    geom_point(data=f0.sum[collection.date!="2019-01-25"&food=="fruit"], aes(x=collection.date, y=prop.d7, color="fruit-reared\noutdoor F0s\nstage 8"))+
    geom_line(data=f0.sum[collection.date!="2019-01-25"&food=="fruit"], aes(x=collection.date, y=prop.d7, color="fruit-reared\noutdoor F0s\nstage 8"))+
    labs(x=NULL, y="Proportion of flies with\ndiapause-like ovaries", color="", linetype="") +
    geom_errorbar(data=f0.sum[collection.date!="2019-01-25"&food=="fruit"], aes(x=collection.date, ymin=prop.d9-se.prop.d9, ymax=prop.d9+se.prop.d9, color="fruit-reared\noutdoor F0s\nstage 10"), width=0.1)+
    geom_errorbar(data=f0.sum[collection.date!="2019-01-25"&food=="fruit"], aes(x=collection.date, ymin=prop.d7-se.prop.d7, ymax=prop.d7+se.prop.d7, color="fruit-reared\noutdoor F0s\nstage 8"), width=0.1)+
    geom_point(data=f0.sum[food=="cornmeal-mollasses"], aes(x=collection.date, y=prop.d9, color = "molasses-reared\noutdoor F0s\nstage 10"))+
    geom_point(data=f0.sum[food=="cornmeal-mollasses"], aes(x=collection.date, y=prop.d7, color = "molasses-reared\noutdoor F0s\nstage 8"))+
        
    geom_errorbar(data=f0.sum[food=="cornmeal-mollasses"], aes(x=collection.date, ymin=prop.d9-se.prop.d9, ymax=prop.d9+se.prop.d9, color = "molasses-reared\noutdoor F0s\nstage 10"), width=0.1)+
        geom_errorbar(data=f0.sum[food=="cornmeal-mollasses"], aes(x=collection.date, ymin=prop.d7-se.prop.d7, ymax=prop.d7+se.prop.d7, color = "molasses-reared\noutdoor F0s\nstage 8"), width=0.1)+
    #geom_hline(aes(yintercept=dat.sum.date[collection.date=="2019-01-25", prop.d9], color="lab-reared indoor F0s", linetype="dashed"))+
    geom_hline(aes(yintercept = f0.sum[collection.date=="2019-01-25", prop.d9], linetype="indoor-reared F0s"), linetype="dashed", color="grey50")+
    geom_hline(aes(yintercept = f0.sum[collection.date=="2019-01-25", prop.d7], linetype="indoor-reared F0s"), linetype="dashed", color="grey80")+
    
    scale_color_manual(values=c("forestgreen", "darkolivegreen3", "saddlebrown", "peru"))+
    scale_x_date(date_breaks = "1 month", date_labels="%b", limits = as.Date(c('2018-06-15','2018-12-15')))

#stats for figure legend
summary(glmer(diapause.bin~as.factor(collection.date)+(1|Cage), data=f0[collection.date!="2019-01-25"&collection.date!="2018-10-15"], family='binomial'))

summary(glmer(diapause.bin~as.factor(collection.date)+(1|Cage), data=f0[collection.date=="2018-10-16"|collection.date=="2018-10-15"], family='binomial'))

summary(glmer(diapause.bin9~as.factor(collection.date)+(1|Cage), data=f0[collection.date!="2019-01-25"&collection.date!="2018-10-15"], family='binomial'))

summary(glmer(diapause.bin9~as.factor(collection.date)+(1|Cage), data=f0[collection.date=="2018-10-16"|collection.date=="2018-10-15"], family='binomial'))


f2<-fread("~/Box Sync/manuscripts/diapause gwas/dryad_upload/mesocosm_F2.csv") #on dryad
f2[eggs==0&ovary<8, diapause.bin:=1]
f2[eggs>0|(eggs==0&ovary>=8), diapause.bin:=0]
f2[,date:=as.Date(date)]

f2.sum=f2[,.(prop.d9=sum(diapause.bin9, na.rm=T)/.N,
             prop.d7=sum(diapause.bin, na.rm=T)/.N,
             med.eggs=as.double(median(as.integer(eggs), na.rm=T)),
             mad.eggs=mad(x=eggs, constant=1/quantile(eggs,.75, na.rm=T)),
             n=.N), 
          .(box,location,date, temp)]
f2.sum[,se.prop.d9:=sqrt(prop.d9*(1-prop.d9)/n)]
f2.sum[,se.prop.d7:=sqrt(prop.d7*(1-prop.d7)/n)]


c<-ggplot()+ geom_point(data=f2.sum, aes(x=date, y=prop.d9, color=location, group=interaction(location,temp)))+
    geom_line(data=f2.sum, aes(x=date, y=prop.d9, color=location, linetype=temp, group=interaction(location, temp)))+
    geom_errorbar(data=f2.sum, aes(x=date, ymin=prop.d9-se.prop.d9, ymax=prop.d9+se.prop.d9, color = location), width=0.1)+ labs(x=NULL, y="Proportion of flies\nin diapause", color="Location", linetype="Temperature") +lims(y=c(0,1))+
    scale_color_manual(values=c("grey50","forestgreen" ))+
    scale_x_date(date_breaks = "1 month", date_labels="%b", limits = as.Date(c('2018-06-15','2018-12-15')))

b<-ggplot()+ geom_point(data=f2.sum, aes(x=date, y=prop.d7, color=location, group=interaction(location,temp)))+
    geom_line(data=f2.sum, aes(x=date, y=prop.d7, color=location, linetype=temp, group=interaction(location, temp)))+
    geom_errorbar(data=f2.sum, aes(x=date, ymin=prop.d7-se.prop.d7, ymax=prop.d7+se.prop.d7, color = location), width=0.1)+ labs(x=NULL, y="Proportion of flies\nin diapause", color="Location", linetype="Temperature") +lims(y=c(0,.5))+
    scale_color_manual(values=c("grey80","darkolivegreen3" ))+
    scale_x_date(date_breaks = "1 month", date_labels="%b", limits = as.Date(c('2018-06-15','2018-12-15')))



library(lme4)

summary(glmer(diapause.bin~as.factor(date)+(1|cage), data=f2[location=="indoor"&box==46], family='binomial'))
summary(glmer(diapause.bin~as.factor(date)+(1|cage), data=f2[location=="outdoor"&box==46], family='binomial'))
summary(glmer(diapause.bin~as.factor(date)+(1|cage), data=f2[location=="indoor"&box==15], family='binomial'))
summary(glmer(diapause.bin~as.factor(date)+(1|cage), data=f2[location=="outdoor"&box==15], family='binomial'))

summary(glmer(diapause.bin9~as.factor(date)+(1|cage), data=f2[location=="indoor"&box==46], family='binomial'))
summary(glmer(diapause.bin9~as.factor(date)+(1|cage), data=f2[location=="outdoor"&box==46], family='binomial'))
summary(glmer(diapause.bin9~as.factor(date)+(1|cage), data=f2[location=="indoor"&box==15], family='binomial'))
summary(glmer(diapause.bin9~as.factor(date)+(1|cage), data=f2[location=="outdoor"&box==15], family='binomial'))




w<-fread("~/Dropbox/Helen-Priscilla/weather_data.csv") #on dryad
w[,date2:=mdy(date)]
w[,jdate:=yday(date2)]

w.melt<-melt(w, id.vars=c("date", "date2", "jdate"), measure.vars = c("High_temp", "Avg_temp", "Low_temp"), variable.name = "temp_type", value.name = "temperature")

w.melt[,tempC:=(temperature-32)*(5/9)]
w.melt[temp_type=="High_temp", temp_type:="High"]
w.melt[temp_type=="Avg_temp", temp_type:="Average                              "]
w.melt[temp_type=="Low_temp", temp_type:="Low"]


tempplot<-ggplot(data=w.melt, aes(x=date2, y=tempC, color=temp_type, group=temp_type))+
    geom_line()+
    scale_x_date(date_breaks = "1 month", date_labels="%b", limits = as.Date(c('2018-06-15','2018-12-15')))+
    labs(x="Date", y="\nTemperature °C\n", color="")+
    geom_hline(yintercept=0, color="gray")


pdf("~/Desktop/Figure_6.pdf", height=8, width=8)

plot_grid(a, b, c, tempplot, align = "v", nrow = 4, rel_heights = c(.3, .3, .3, .2), labels=c("A", "B", "C", "D"))

dev.off()

######################################
###### FIGURE 7 : Africa #############
######################################

load("/scratch/pae3g/revisions/all_lasso_for_alan.Rdata") #on dryad

zi<-fread("/scratch/pae3g/revisions/evolution/fst_ZI_AT_gr_12_fall.txt") #on dryad

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
lasso.zi.sum[,top:="lasso"]


lasso.zi.sum.melt<-melt(lasso.zi.sum, id.vars=c("pop", "pheno", "GRM", "perm" ,"top"), measure.vars=c("med", "prop.01", "prop.05", "prop.1", "prop.2"))
setnames(lasso.zi.sum.melt, "GRM", "draw")


load("/scratch/pae3g/revisions/evolution/zambia_universal_threshold.Rdata") # on dryad
zambia.melt<-melt(zambia, id.vars=c("pop", "pheno", "draw", "perm", "top"), measure.vars=c("med", "prop.01", "prop.05", "prop.1", "prop.2"))

zi.all<-rbind(zambia.melt, lasso.zi.sum.melt)

zi.all[perm!=0&pop=="both", group:="both-permuted"]
zi.all[perm==0&pop=="both", group:="both-observed"]
zi.all[perm==0&pop=="A", group:="A-observed"]
zi.all[perm==0&pop=="B", group:="B-observed"]
zi.all[perm!=0&pop=="A", group:="A-permuted"]
zi.all[perm!=0&pop=="B", group:="B-permuted"]

zi.all[,pheno2:=ifelse(pheno=="diapause.bin", "st. 8", "st. 10")]
zi.all[,pheno2:=factor(zi.all$pheno2, levels=c("st. 8", "st. 10"))]


zi.all.sum<-zi.all[,.(med=median(value, na.rm=T), q.025=quantile(value, 0.025, na.rm=T), q.975=quantile(value, .975, na.rm=T)), .(group, pheno, pheno2, pop,top, variable)]
zi.all.sum[,pheno2:=factor(zi.all.sum$pheno2, levels=c("st. 8", "st. 10"))]


zi.all.sum[top==-5, th:="Top 0.001%"]
zi.all.sum[top==-4, th:="Top 0.01%"]
zi.all.sum[top==-3, th:="Top 0.1%"]
zi.all.sum[top==-2, th:="Top 1%"]
zi.all.sum[top==-1, th:="Top 10%"]
zi.all.sum[top==0, th:="all SNPs"]
zi.all.sum[top=="lasso", th:="LASSO"]


a.plot<-ggplot(zi.all.sum[variable=="med" & th%in%c("LASSO", "all SNPs", "Top 0.01%", "Top 0.1%")])+
    geom_point(data=zi.all.sum[variable=="med"&th%in%c("LASSO", "all SNPs", "Top 0.01%", "Top 0.1%")], aes(x=pop, y=med, color=group), position=position_dodge(width=0.5))+
    geom_errorbar(data=zi.all.sum[variable=="med"&th%in%c("LASSO", "all SNPs", "Top 0.01%", "Top 0.1%")], aes(x=pop, ymax=q.975, ymin=q.025, color=group), width=0.2, position=position_dodge(width=0.5))+
    labs(x="", y="median frequency of pro-diapause allele in Africa", color="")+
    scale_color_manual(values=c("dodgerblue2", "grey80", "darkorchid", "grey80", "lightseagreen", "grey80"))+
    facet_grid(th~pheno2, scales ="free_y")+
    theme(legend.position = "none")

pdf("/scratch/pae3g/revisions/figures/zambian_main.pdf", height=6, width=5)
a.plot
dev.off()


#############################
###### FIGURE 8 #############
#############################

load("/scratch/pae3g/revisions/evolution/ihs_dgrp_universal_lasso.Rdata") # on dryad

ihs.dgrp<-copy(p)
load("/scratch/pae3g/revisions/evolution/ihs_northern_universal_lasso.Rdata") # on dryad

ihs.north<-copy(p)

a<-rbind(ihs.dgrp, ihs.north)
lasso.ihs.sum<-a[,.(min.ihs=min(IHS, na.rm=T),
                    med.ihs=median(IHS, na.rm=T), 
                    max.ihs=max(IHS, na.rm=T),
                    n=.N), .(pop, pheno, draw, perm, test, model)]

lasso.ihs.sum[,top:="lasso"]
lasso.ihs.sum[test=="dgrp", test:="DGRP"]
lasso.ihs.sum[test=="north", test:="Northern"]



load("/scratch/pae3g/revisions/evolution/ihs_dgrp_universal_threshold.Rdat") # on dryad
load("/scratch/pae3g/revisions/evolution/ihs_north_universal_threshold.Rdat") # on dryad

j[,test:="DGRP"]
k[,test:="Northern"]

i<-rbind(j,k)

ihs<-rbind(lasso.ihs.sum, i)
ihs.melt<-melt(ihs, measure.vars=c("min.ihs", "med.ihs", "max.ihs"), id.vars=c("pop", 'pheno', "draw", "perm", "test", "top"))

ihs.melt[perm!=0&pop=="both", group:="both-permuted"]
ihs.melt[perm==0&pop=="both", group:="both-observed"]
ihs.melt[perm==0&pop=="A", group:="A-observed"]
ihs.melt[perm==0&pop=="B", group:="B-observed"]
ihs.melt[perm!=0&pop=="A", group:="A-permuted"]
ihs.melt[perm!=0&pop=="B", group:="B-permuted"]

ihs.melt[,pheno2:=ifelse(pheno=="diapause.bin", "st. 8", "st. 10")]
ihs.melt[,pheno2:=factor(ihs.melt$pheno2, levels=c("st. 8", "st. 10"))]


ihs.melt.sum<-ihs.melt[,.(med=median(value, na.rm=T), q.025=quantile(value, 0.025, na.rm=T), q.975=quantile(value, .975, na.rm=T)), .(group, pheno, pheno2, pop,top, variable, test)]
ihs.melt.sum[,pheno2:=factor(ihs.melt.sum$pheno2, levels=c("st. 8", "st. 10"))]


ihs.melt.sum[top==-5, th:="Top 0.001%"]
ihs.melt.sum[top==-4, th:="Top 0.01%"]
ihs.melt.sum[top==-3, th:="Top 0.1%"]
ihs.melt.sum[top==-2, th:="Top 1%"]
ihs.melt.sum[top==-1, th:="Top 10%"]
ihs.melt.sum[top==0, th:="all SNPs"]
ihs.melt.sum[top=="lasso", th:="LASSO"]

ihs.melt.sum[,facet:=paste(pheno2, test, sep=": ")]
ihs.melt.sum[,facet:=factor(ihs.melt.sum$facet, levels=c("st. 8: DGRP", "st. 10: DGRP", "st. 8: Northern", "st. 10: Northern"))]


a.plot<-ggplot(ihs.melt.sum[test=="DGRP"& variable=="med.ihs" & th%in%c("LASSO", "all SNPs", "Top 0.01%", "Top 0.1%")])+
    geom_point(data=ihs.melt.sum[test=="DGRP"& variable=="med.ihs"&th%in%c("LASSO", "all SNPs", "Top 0.01%", "Top 0.1%")], aes(x=pop, y=med, color=group), position=position_dodge(width=0.5))+
    geom_errorbar(data=ihs.melt.sum[test=="DGRP"& variable=="med.ihs"&th%in%c("LASSO", "all SNPs", "Top 0.01%", "Top 0.1%")], aes(x=pop, ymax=q.975, ymin=q.025, color=group), width=0.2, position=position_dodge(width=0.5))+
    labs(x="", y="median iHS", color="", title="DGRP")+
    scale_color_manual(values=c("dodgerblue2", "grey80", "darkorchid", "grey80", "lightseagreen", "grey80"))+
    facet_grid(th~pheno2, scales ="free_y")+
    theme(legend.position = "none")


b.plot<-ggplot(ihs.melt.sum[test=="Northern"& variable=="med.ihs" & th%in%c("LASSO", "all SNPs", "Top 0.01%", "Top 0.1%")])+
    geom_point(data=ihs.melt.sum[test=="Northern"&variable=="med.ihs"&th%in%c("LASSO", "all SNPs", "Top 0.01%", "Top 0.1%")], aes(x=pop, y=med, color=group), position=position_dodge(width=0.5))+
    geom_errorbar(data=ihs.melt.sum[test=="Northern"&variable=="med.ihs"&th%in%c("LASSO", "all SNPs", "Top 0.01%", "Top 0.1%")], aes(x=pop, ymax=q.975, ymin=q.025, color=group), width=0.2, position=position_dodge(width=0.5))+
    labs(x="", y=" ", color="", title="Northern")+
    scale_color_manual(values=c("dodgerblue2", "grey80", "darkorchid", "grey80", "lightseagreen", "grey80"))+
    facet_grid(th~pheno2, scales ="free_y")+
    theme(legend.position = "none")
pdf("/scratch/pae3g/revisions/figures/ihs_main.pdf", height=8, width=8)
plot_grid(a.plot, b.plot, nrow=1, labels=c("A", "B"))
dev.off()

####################################################
########## Fig S1: PCA in parents ##################
####################################################


pca<-fread("/nv/vol186/bergland-lab/Priscilla/parents_PCA.txt") #on dryad
pca.season<-fread("/nv/vol186/bergland-lab/Priscilla/parents_season_PCA.txt") #on dryad
pca.lat<-fread("/nv/vol186/bergland-lab/Priscilla/parents_northsouth_PCA.txt") #on dryad

pca[,chrs:=factor(pca$chrs, levels=c("all", "2L", "2R", "3L", "3R", "X"))]
pca.season[,chrs:=factor(pca.season$chrs, levels=c("all", "2L", "2R", "3L", "3R", "X"))]
pca.lat[,chrs:=factor(pca.lat$chrs, levels=c("all", "2L", "2R", "3L", "3R", "X"))]

pca.lat[geography=="Maine", loc:="North"]
pca.lat[geography=="Penn-spring", loc:="North"]
pca.lat[geography=="Penn-fall", loc:="North"]
pca.lat[geography=="New York", loc:="North"]
pca.lat[geography=="Southeast", loc:="South"]
pca.lat[geography=="Carribean", loc:="South"]


#with all chromosomes

a<-ggplot(pca, aes(x=pc1, y=pc2, color=geography))+
    geom_point()+
    labs(x="PC 1", y="PC 2", color="founding line")+
    facet_grid(.~chrs)

b<-ggplot(pca.lat, aes(x=pc1, y=pc2, color=loc))+
    geom_point()+
    labs(x="PC 1", y="PC 2", color="founding line")+
    facet_grid(.~chrs)

c<-ggplot(pca.season[(geography=="Penn-spring"|geography=="Penn-fall")], aes(x=pc1, y=pc2, color=geography))+
    geom_point()+
    labs(x="PC 1", y="PC 2", color="founding line")+
    facet_grid(.~chrs)

d<-ggplot()

pdf("/scratch/pae3g/revisions/figures/parent_line_PCA_allchrs.pdf", height=6, width=10)
plot_grid(a,b,c, nrow=3, labels = c("A", "B", "C"), align="vh",  axis="lrbt")
dev.off()


#################################################
####### Fig S2: Hybrid and parent PCA ###########
#################################################


load("/scratch/pae3g/revisions/hybrid_parent_PCA.Rdata")#on dryad
load("/scratch/pae3g/revisions/parents_hybrids_karyotype_calls.Rdata") #on dryad
kary.ag.melt=melt(kary.ag, measure.vars=c("chr2L", "chr2R", "chr3L", "chr3R"), id.vars="sample.id")

pcaOut=merge(pcaOut2, kary.ag.melt, by=c("sample.id", "chrs"), all=TRUE)
pcaOut[is.na(kary), kary:="std"]
pcaOut[chrs=="chr2L", chr:="2L"]
pcaOut[chrs=="chr3L", chr:="3L"]
pcaOut[chrs=="chr2R", chr:="2R"]
pcaOut[chrs=="chr3R", chr:="3R"]
pcaOut[chrs=="chrX", chr:="X"]
pcaOut[chrs=="all", chr:="all"]
names(kary.ag.melt)=c("sample.id", "chrs", "kary")

pcaOut[,kary2:=ifelse(kary=="std", "Standard", "Inversion")]
pcaOut[,chr:=factor(chr, levels=c("all", "2L", "2R", "3L", "3R", "X"))]


#reorder karyotype levels
pcaOut[,kary:=factor(kary, levels=c("std", "std;In(2L)t" , "std;In(2R)Ns" , "std;In(3R)C", "std;In(3R)Mo" , "std;In(3R)Payne","In(2L)t","In(2R)Ns","In(3R)C", "In(3R)Mo","In(3R)Payne", "In(3R)Payne;In(3R)C", "In(3R)Payne;In(3R)Mo" ))]



pdf("/scratch/pae3g/revisions/figures/pca_kary_hybrids_parents.pdf", height=6, width=10)
ggplot(pcaOut, aes(x=pc1, y=pc2, color=group, shape=kary2, alpha=type))+
    geom_point()+
    facet_wrap(~chr, scales="free")+
    labs(color="group", shape="karyotype", alpha="population")+
    scale_color_manual(values=c( "steelblue2", "mediumpurple2","navyblue","purple4"  ))+
    scale_alpha_manual(values=c(0.2, 1))
dev.off()



###################################
#### FIGURE S3 ####################
###################################

#average calibrated temperature of each box, color coded by photoperiod

library(data.table)
library(foreach)
library(cowplot)
library(lubridate)
library(ggbeeswarm)
library(foreach)

env<-fread("~/Box Sync/hybridSwarm/all_mapping_environmental_data.txt")

#pull out time info
env[,time:=ymd_hms(TimeStamp)]
env[,hour:=hour(time)]
env[,min:=minute(time)]
env[,day:=yday(time)]
env[,adj.day:=day-min(env$day, na.rm=T)]

#merge with info about boxes and temperature calibration
env<-merge(env, info, by="name")
env<-merge(env, cal, by="box")

env[photoperiod==9, pp:="9L:15D"]
env[photoperiod==11, pp:="11L:13D"]
env[photoperiod==13, pp:="13L:11D"]
env[photoperiod==15, pp:="15L:9D"]
env[,pp:=factor(env$pp, levels=c("9L:15D", "11L:13D", "13L:11D", "15L:9D"))]


#calculate calibrated temperature
env[,temp.rack.cal:=MCP9808Temp-mcp.rack.diff]

#summarize the data to get an average for each hour of the day
env.hoursum<-env[,.(mean.temp=mean(temp.rack.cal), sd.temp=sd(temp.rack.cal)), .(hour, box_num, pp)]
env.hoursum[,ci95:=1.96*sd.temp]

#plot hourly mean temperatures for every box, facet by photoperiod

a<-ggplot(env.hoursum, 
          aes(x=hour, 
              y=mean.temp, 
              fill=box_num))+
    facet_grid(pp~.)+ 
    theme(legend.position="none")+
    geom_line()+
    labs(x="Hour of Day", 
         y="Mean temperature °C")


#subset data to one average recording per hour
env.sum<-env[,.(mean.temp=mean(temp.rack.cal), sd.temp=sd(temp.rack.cal)), .(hour, box_num, pp, day, adj.day)]
env.sum[,ci95:=1.96*sd.temp]
env.sum[,avg.time:=adj.day+hour/24]

#ggplot(env.sum, aes(x=avg.time, y=mean.temp, color=as.factor(box_num), group(box_num)))+geom_line()+ theme(legend.position="none")
#summarize temperatures by lights on and lights off

env<-env[Lights=="Off"|Lights=="on, main"]
lights.sum<-env[,.(mean.temp=mean(temp.rack.cal), sd.temp=sd(temp.rack.cal)), .(box, Lights, pp)]


#plot temperature of each box with lights on and off

b<-ggplot(lights.sum, 
          aes(x=Lights, 
              y=mean.temp, 
              group=box, 
              color=as.factor(pp)))+
    geom_point()+geom_line()+ 
    scale_x_discrete(labels=c("Lights off", "Lights on"))+
    labs(x=NULL, 
         y="Mean temperature °C", 
         color="Photoperiod")#+geom_errorbar(aes(ymin=mean.temp-sd.temp, ymax=mean.temp+sd.temp, width=0.05))

#plot histogram of delta temps

lights.wide<-dcast(lights.sum, box+pp~Lights, value.var="mean.temp")
setnames(lights.wide, "on, main", "On")

lights.wide[,delta:=On-Off]

c<-ggplot(lights.wide, 
          aes(x=as.factor(pp), 
              y=delta, 
              color=pp))+
    geom_boxplot()+
    geom_beeswarm()+
    labs(x="Photoperiod", 
         y="On - off (°C)")+
    theme(axis.text.x=element_text(angle=45,hjust=1))

#final temperature plots
right<-plot_grid(b+theme(legend.position="none"), 
                 c+theme(legend.position="none"),  
                 rel_heights = c(0.6, 0.4), 
                 labels=c("B", "C"), 
                 nrow=2)
legend <- get_legend(b)

pdf("~/Box Sync/manuscripts/FigureS1.pdf", height=8, width=8)
plot_grid(a, right, legend, 
          rel_widths = c(0.4, 0.4, 0.2), 
          labels=c("A", ""), 
          nrow=1)
dev.off()



#################################
#### FIGURE S4 ##################
#################################


#supplemental figure showing generation and swarm with each phenotype continues from the data table generated in Figure 2

#melt phenotypes into long format
p.melt<-melt(p, id.vars=c("id", "temp.rack.cal", "photoperiod", "generation", "swarm"), measure.vars=c("diapause.bin9", "diapause.bin", "eggP"), variable.name="pheno", value.name="diapause")
p.melt[,pheno:=factor(p.melt$pheno, levels=c("eggP", "diapause.bin9", "diapause.bin"))]

p.melt[pheno=="eggP", stage:="Stage 14"]
p.melt[pheno=="diapause.bin", stage:="Stage 8"]
p.melt[pheno=="diapause.bin9", stage:="Stage 10"]

p.melt[,stage:=factor(p.melt$stage, levels=c("Stage 8", "Stage 10", "Stage 14"))]

p.melt[generation==4,diapause.gen:=as.numeric(diapause)+.01]
p.melt[generation==5,diapause.gen:=as.numeric(diapause)-.01]

p.melt[swarm=="A",diapause.swarm:=as.numeric(diapause)+.01]
p.melt[swarm=="B",diapause.swarm:=as.numeric(diapause)-.01]


#Generation plot
a.1<-ggplot(p.melt, aes(x=temp.rack.cal, y=diapause, color=as.factor(generation)))+
    binomial_smooth()+
    geom_point(data=p.melt, aes(x=temp.rack.cal, y=diapause.gen,color=as.factor(generation)), alpha=0.1)+
    labs(color="Generation", x="Temperature °C", y="Probability of diapause")+facet_grid(.~stage)

#Cage plot
b.1<-ggplot(p.melt, aes(x=temp.rack.cal, y=diapause, color=as.factor(swarm)))+
    binomial_smooth()+
    geom_point(data=p.melt, aes(x=temp.rack.cal, y=diapause.swarm, color=as.factor(swarm)), alpha=0.1)+
    labs(color="Population", x="Temperature °C", y="Probability of diapause")+facet_grid(.~stage)
#Final plot: effect of generation and cage on phenotypes

pdf("~/Box Sync/manuscripts/diapause gwas/R1/figureS2.pdf", height=6, width=8)
plot_grid(a.1, b.1, nrow=2, labels=c("A", "B"))
dev.off()




#################################
#### FIGURE S5 Yeast ############
#################################

library(data.table)
library(cowplot)


data<-fread("~/Box Sync/manuscripts/diapause gwas/dryad_upload/yeast_supplementation.csv")

p.melt<-melt(data, id.vars=c( "temp.rack.cal", "box", "treatment", "cage"), measure.vars=c("diapause9", "diapause", "eggp"), variable.name="pheno", value.name="diapause")

p.melt[pheno=="eggp", stage:="Stage 14"]
p.melt[pheno=="diapause", stage:="Stage 8"]
p.melt[pheno=="diapause9", stage:="Stage 10"]

p.melt[,stage:=factor(p.melt$stage, levels=c("Stage 8", "Stage 10", "Stage 14"))]
p.melt[treatment=="yeast", diapause.adjust:=as.numeric(diapause)+.02]
p.melt[treatment=="no yeast", diapause.adjust:=as.numeric(diapause)-.02]


pdf("~/Box Sync/manuscripts/diapause gwas/R1/FigureS3.pdf", heigh=3.5, width=8)
ggplot(p.melt, aes(x=temp.rack.cal,y=diapause, color=treatment, linetype=cage))+
        binomial_smooth()+
        labs(x="Temperature (°C)",y="Probability of diapause")+
        scale_y_continuous(breaks=c(0,1))+
        geom_jitter(data=p.melt, aes(x=temp.rack.cal, y=diapause.adjust, color=treatment, alpha=0.05))+
    labs(color=NULL, alpha=NULL)+
    facet_wrap(~stage)
dev.off()

summary(aov(glm(diapause~temp.rack.cal+cage+yeast, family=binomial, data=data)))
summary(aov(glm(diapause9~temp.rack.cal+cage+yeast, family=binomial, data=data)))
summary(aov(glm(eggp~temp.rack.cal+cage+yeast, family=binomial, data=data)))




###############################################
## Figure S6: chromosome painting  ############
###############################################

keep.paths<-fread("/scratch/pae3g/evolution/all_reconstructed_haplotypes_melted.txt")
keep.paths[path.length<1000000, line:="UNKNOWN"]

#assign arbitrary y-values for plot
keep.paths<-keep.paths[order(swarm)]
keep.paths[,index:=rleid(sample.id)]
keep.paths[haplotype=="par1", y1:=as.numeric(index)]
keep.paths[haplotype=="par1", y2:=index+.4]
keep.paths[haplotype=="par2", y1:=index+.4]
keep.paths[haplotype=="par2", y2:=index+.8]


# FINAL PLOT

pdf("/scratch/pae3g/reconstructions.pdf", height=8.75, width=7.5)
ggplot()+geom_rect(data=keep.paths[line!="UNKNOWN"], 
                   mapping=aes(xmin=cons.start/1000000, 
                               xmax=cons.stop/1000000, 
                               ymin=y1, 
                               ymax=y2, 
                               fill=line))+ 
    theme(legend.position="bottom",  
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(), 
          axis.line.y=element_blank(),
          legend.text =element_text(size=8),
          legend.key.size = unit(.1, "in"))+
    labs(x="Position, Mb") +
    facet_grid(swarm~chromosome, scales="free") +
    geom_rect(data=keep.paths[line=="UNKNOWN"], mapping=aes(xmin=cons.start/1000000, xmax=cons.stop/1000000, ymin=y1, ymax=y2))
dev.off()





########################################
### FIGURE S7  GIF by chr ##############
########################################


files<-fread("/scratch/pae3g/genome-reconstruction/universal_input.txt")

p<-foreach(pop=files$V1[1:3000], phenotype=files$V2[1:3000],draw=files$V3[1:3000], perm=files$V4[1:3000], model=files$V6[1:3000], .errorhandling="remove")%dopar%{
    load(paste0("/scratch/pae3g/revisions/lasso/", "GIFbychr_", phenotype, "_draw",draw, "_perm", perm, "_", model, "_pop", pop, ".Rdata"))
    return(x)
}

p<-rbindlist(p)

save (p, file= "/scratch/pae3g/revisions/GIFbychr.Rdata") # on dryad

p[perm!=0&pop=="both", group:="both-permuted"]
p[perm==0&pop=="both", group:="both-observed"]
p[perm==0&pop=="A", group:="A-observed"]
p[perm==0&pop=="B", group:="B-observed"]
p[perm!=0&pop=="A", group:="A-permuted"]
p[perm!=0&pop=="B", group:="B-permuted"]
p[,pheno2:=ifelse(pheno=="diapause.bin", "st. 8", "st. 10")]
p[,pheno2:=factor(p$pheno2, levels=c("st. 8", "st. 10"))]


p.sum<-p[,.(med=median(lambda.05), q.025=quantile(lambda.05, 0.025), q.975=quantile(lambda.05, .975)), .(group, pheno, pheno2, chr, pop)]
p.sum[,pheno2:=factor(p.sum$pheno2, levels=c("st. 8", "st. 10"))]

a.plot<-ggplot(p.sum)+
    geom_point(data=p.sum, aes(x=pop, y=med, color=group), position=position_dodge(width=0.5))+
    geom_errorbar(data=p.sum, aes(x=pop, ymax=q.975, ymin=q.025, color=group), width=0.2, position=position_dodge(width=0.5))+
    labs(x="", y="GIF", color="")+
    scale_color_manual(values=c("dodgerblue2", "grey80", "darkorchid", "grey80", "lightseagreen", "grey80"))+
    theme(legend.position = "none")+
    facet_grid(pheno2~chr, scales ="free_y")


pdf("/scratch/pae3g/revisions/figures/GIF_by_chr.pdf", height=5, width=8)
a.plot
dev.off()


##########################
#### FIGURE S8 ##########
##########################


#make LD plots of LASSO SNPs for 10 imputations and 10 random permutations
library(data.table)
library(foreach)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(viridis)
library(SNPRelate)
library(doMC)
registerDoMC(10)


files<-fread("/scratch/pae3g/genome-reconstruction/universal_input.txt")
files<-files[V6=="nonloco"&V1=="both"&V2=="diapause.bin9"&((V4>=1000&V4<1010)| (V4==0&V3<=10))]

#read in files
ld<-foreach(pop=files$V1, phenotype=files$V2,draw=files$V3, perm=files$V4, model=files$V6, .errorhandling="remove")%dopar%{
    print(paste(draw, perm, sep=","))
    load(paste0("/scratch/pae3g/revisions/lasso/", "lasso_", phenotype, "_draw",draw, "_perm", perm, "_", model, "_pop", pop, ".lassoSites.Rdata"))
    sites[, chr:=tstrsplit(site,split="_")[[1]]]
    sites[, pos:=as.numeric(tstrsplit(site,split="_")[[2]])]
    geno<-snpgdsOpen(paste0("/scratch/pae3g/genome-reconstruction/final2_draw", draw, "_replaced.vcf.gds"))
    a<-snpgdsSNPList(gdsobj=geno)
    info<-data.table(snp.id=a$snp.id,
                     chr=a$chromosome,
                     pos=a$pos)
    sites<-merge(sites, info, by=c("chr", "pos"))
    ldmat<-snpgdsLDMat(geno, slide=-1,verbose=F,snp.id =sites$snp.id , method='composite')$LD
    m<-data.table(melt(ldmat))
    m[,perm:=perm]
    m[,draw:=draw]
    snpgdsClose(geno)
    return(m)
}

ld<-rbindlist(ld)
ld[perm==0, type:="observed"]
ld[perm!=0, type:="permuted"]
ld[perm>=1000, draw:=perm-999]


pdf("/scratch/pae3g/revisions/figures/ld_heatmaps_lasso.pdf", height=7, width=8)
a<-ggplot(ld[type=="observed"], aes(x=Var1, y=Var2, fill=value*value))+
    geom_tile()+
    scale_fill_viridis(option="viridis")+ 
    theme(line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())+
    labs( fill=expression("R"^2))+
    facet_wrap(~draw, scales="free", nrow = 2)
b<-ggplot(ld[type=="permuted"], aes(x=Var1, y=Var2, fill=value*value))+
    geom_tile()+
    scale_fill_viridis(option="viridis")+ 
    theme(line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())+
    labs( fill=expression("R"^2))+
    facet_wrap(~draw, scales="free", nrow = 2)
plot_grid(a,b,nrow=2, labels=c("A", "B"))
dev.off()

###########################################
### FIGURE S9  ############################
###########################################


files<-fread("/scratch/pae3g/genome-reconstruction/universal_input.txt")

p<-foreach(pop=files$V1[1:3000], phenotype=files$V2[1:3000],draw=files$V3[1:3000], perm=files$V4[1:3000], model=files$V6[1:3000], .errorhandling="remove")%dopar%{
    load(paste0("/scratch/pae3g/revisions/lasso/", "lasso_", phenotype, "_draw",draw, "_perm", perm, "_", model, "_pop", pop, ".lassoSites.Rdata"))
    return(sites)
}

p<-rbindlist(p)

p[, chr:=tstrsplit(site,split="_")[[1]]]
p[, pos:=as.numeric(tstrsplit(site,split="_")[[2]])]

p<-p[chr%in%(c("2L", "2R", "3L", "3R", "X"))]

p.sum<-p[,.(n=.N), .(pop, pheno, loco, GRM, perm)]

p.sum[perm!=0&pop=="both", group:="both-permuted"]
p.sum[perm==0&pop=="both", group:="both-observed"]
p.sum[perm==0&pop=="A", group:="A-observed"]
p.sum[perm!=0&pop=="A", group:="A-permuted"]

p.sum[perm==0&pop=="B", group:="B-observed"]
p.sum[perm!=0&pop=="B", group:="B-permuted"]

p.sum[,pheno2:=ifelse(pheno=="diapause.bin", "st. 8", "st. 10")]
p.sum[,pheno2:=factor(p.sum$pheno2, levels=c("st. 8", "st. 10"))]


p.sum.sum<-p.sum[,.(med=median(n), q.025=quantile(n, 0.025), q.975=quantile(n, .975)), .(group, pheno, pheno2, pop)]
p.sum.sum[,pheno2:=factor(p.sum.sum$pheno2, levels=c("st. 8", "st. 10"))]

a.plot<-ggplot(p.sum)+
    geom_point(data=p.sum.sum, aes(x=pop, y=med, color=group), position=position_dodge(width=0.5))+
    geom_errorbar(data=p.sum.sum, aes(x=pop, ymax=q.975, ymin=q.025, color=group), width=0.1, position=position_dodge(width=0.5))+
    labs(x="", y="Number of LASSO SNPs", color="")+
    scale_color_manual(values=c("dodgerblue2", "grey80", "darkorchid", "grey80", "lightseagreen", "grey80"))+
    facet_grid(.~pheno2, scales ="free_y")+
    theme(legend.position = "none")

pdf("/scratch/pae3g/revisions/figures/number_of_lasso_snps.pdf", height=5, width=6)
a.plot
dev.off()



################################
##### FIGURE S10 ###############
################################

 #see "phenotype_population_overlaps.R" script for analysis and plotting of this figure


#####################################
## FIGURE S11: Cpo and Tim ##########
#####################################

library(data.table)
library(SNPRelate)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

tim.geno<-fread("/scratch/pae3g/tim_genotypes.csv") #on dryad

phenos<-fread("/nv/vol186/bergland-lab/Priscilla/phenos_wolbachia.txt") #on dryad
phenos[photoperiod==9, pp:="9L:15D"]
phenos[photoperiod==11, pp:="11L:13D"]
phenos[photoperiod==13, pp:="13L:11D"]
phenos[photoperiod==15, pp:="15L:9D"]
phenos[,pp:=factor(phenos$pp, levels=c("9L:15D", "11L:13D", "13L:11D", "15L:9D"))]
tim.genos<-merge(phenos, tim.geno, by="id")


#open up hybrid file to pull cpo genotypes
geno<-snpgdsOpen("/scratch/pae3g/oldscratch_recovered/genome-reconstruction/final3.vcf.gds") #gds of final2.vcf, which is on dryad

snp.dt <- data.table(snp.id=read.gdsn(index.gdsn(geno, "snp.id")),
                     chr=read.gdsn(index.gdsn(geno, "snp.chromosome")),
                     pos=read.gdsn(index.gdsn(geno, "snp.position")),
                     alleles=read.gdsn(index.gdsn(geno, "snp.allele")))

snp.dt[pos==13793588] #snp.id==2313328 #alleles =T/A
y<-snpgdsGetGeno(geno, snp.id="2313328", with.id=T)

cpo<-data.table("cpo.geno"=y$genotype[,1],
                sample.id=y$sample.id)

cpo<-merge(cpo, phenos, by="sample.id")


binomial_smooth <- function(...) {
    geom_smooth(method = "glm", method.args = list(family = "binomial"), se=FALSE)
}

st8.both.tim<-ggplot(data=tim.genos[!is.na(final.geno)], aes(x=temp.rack.cal, y=diapause.bin, color=as.factor(final.geno)))+
    geom_jitter(data=tim.genos[!is.na(final.geno)], aes(x=temp.rack.cal, y=diapause.bin, color=as.factor(final.geno)), height = 0.1, size=0.5)+
    labs(x="", y="Stage 8", color="tim genotype", title="both")+
    binomial_smooth()+
    scale_y_continuous(breaks=c(0,1))+
    theme(legend.position="none")

st8.A.tim<-ggplot(data=tim.genos[!is.na(final.geno)&swarm.y=="A"], aes(x=temp.rack.cal, y=diapause.bin, color=as.factor(final.geno)))+
    geom_jitter(data=tim.genos[!is.na(final.geno)&swarm.y=="A"], aes(x=temp.rack.cal, y=diapause.bin, color=as.factor(final.geno)), height = 0.1, size=0.5)+
    labs(x="", y="", color="tim genotype", title="population A")+
    binomial_smooth()+
    scale_y_continuous(breaks=c(0,1))+
    theme(legend.position="none")

st8.B.tim<-ggplot(data=tim.genos[!is.na(final.geno)&swarm.y=="B"], aes(x=temp.rack.cal, y=diapause.bin, color=as.factor(final.geno)))+
    geom_jitter(data=tim.genos[!is.na(final.geno)&swarm.y=="B"], aes(x=temp.rack.cal, y=diapause.bin, color=as.factor(final.geno)), height = 0.1, size=0.5)+
    labs(x="", y="", color="tim geno", title="population B")+
    binomial_smooth()+
    scale_y_continuous(breaks=c(0,1))+
   scale_color_discrete(labels=c("s/s", "s/ls", "ls/ls"))+
    theme(legend.position=c(0.6, .7))


st10.both.tim<-ggplot(data=tim.genos[!is.na(final.geno)], aes(x=temp.rack.cal, y=diapause.bin9, color=as.factor(final.geno)))+
    geom_jitter(data=tim.genos[!is.na(final.geno)], aes(x=temp.rack.cal, y=diapause.bin9, color=as.factor(final.geno)), height = 0.1, size=0.5)+
    labs(x="", y="Stage 10", color="tim genotype")+
    binomial_smooth()+
    scale_y_continuous(breaks=c(0,1))+
    theme(legend.position="none")

st10.A.tim<-ggplot(data=tim.genos[!is.na(final.geno)&swarm.y=="A"], aes(x=temp.rack.cal, y=diapause.bin9, color=as.factor(final.geno)))+
    geom_jitter(data=tim.genos[!is.na(final.geno)&swarm.y=="A"], aes(x=temp.rack.cal, y=diapause.bin9, color=as.factor(final.geno)), height = 0.1, size=0.5)+
    labs(x="", y="", color="tim genotype")+
    binomial_smooth()+
    scale_y_continuous(breaks=c(0,1))+
    theme(legend.position="none")

st10.B.tim<-ggplot(data=tim.genos[!is.na(final.geno)&swarm.y=="B"], aes(x=temp.rack.cal, y=diapause.bin9, color=as.factor(final.geno)))+
    geom_jitter(data=tim.genos[!is.na(final.geno)&swarm.y=="B"], aes(x=temp.rack.cal, y=diapause.bin9, color=as.factor(final.geno)), height = 0.1, size=0.5)+
    labs(x="", y="", color="tim geno")+
    binomial_smooth()+
    scale_y_continuous(breaks=c(0,1))+
    scale_color_discrete(labels=c("s/s", "s/ls", "ls/ls"))+
    theme(legend.position=c(0.6, .7))
    

st8.both.cpo<-ggplot(data=cpo[!is.na(cpo.geno)], aes(x=temp.rack.cal, y=diapause.bin, color=as.factor(cpo.geno)))+
    geom_jitter(data=cpo[!is.na(cpo.geno)], aes(x=temp.rack.cal, y=diapause.bin, color=as.factor(cpo.geno)), height = 0.1, size=0.5)+
    labs(x="", y="Stage 8", color="cpo genotype")+
    binomial_smooth()+
    scale_y_continuous(breaks=c(0,1))+
    scale_color_discrete(labels=c("A/A", "A/T", "T/T"))+
    theme(legend.position="none")


st8.A.cpo<-ggplot(data=cpo[!is.na(cpo.geno)&swarm=="A"], aes(x=temp.rack.cal, y=diapause.bin, color=as.factor(cpo.geno)))+
    geom_jitter(data=cpo[!is.na(cpo.geno)&swarm=="A"], aes(x=temp.rack.cal, y=diapause.bin, color=as.factor(cpo.geno)), height = 0.1, size=0.5)+
    labs(x="", y="", color="cpo genotype")+
    binomial_smooth()+
    scale_y_continuous(breaks=c(0,1))+
    scale_color_discrete(labels=c("A/A", "A/T", "T/T"))+
    theme(legend.position="none")


st8.B.cpo<-ggplot(data=cpo[!is.na(cpo.geno)&swarm=="B"], aes(x=temp.rack.cal, y=diapause.bin, color=as.factor(cpo.geno)))+
    geom_jitter(data=cpo[!is.na(cpo.geno)&swarm=="B"], aes(x=temp.rack.cal, y=diapause.bin, color=as.factor(cpo.geno)), height = 0.1, size=0.5)+
    labs(x="", y="", color="cpo geno")+
    binomial_smooth()+
    scale_y_continuous(breaks=c(0,1))+
    scale_color_discrete(labels=c("A/A", "A/T", "T/T"))+
    theme(legend.position=c(0.6, 0.7))


st10.both.cpo<-ggplot(data=cpo[!is.na(cpo.geno)], aes(x=temp.rack.cal, y=diapause.bin9, color=as.factor(cpo.geno)))+
    geom_jitter(data=cpo[!is.na(cpo.geno)], aes(x=temp.rack.cal, y=diapause.bin9, color=as.factor(cpo.geno)), height = 0.1, size=0.5)+
    labs(x="Temperature °C", y="Stage 10", color="cpo genotype")+
    binomial_smooth()+
    scale_y_continuous(breaks=c(0,1))+
    scale_color_discrete(labels=c("A/A", "A/T", "T/T"))+
    theme(legend.position="none")


st10.A.cpo<-ggplot(data=cpo[!is.na(cpo.geno)&swarm=="A"], aes(x=temp.rack.cal, y=diapause.bin9, color=as.factor(cpo.geno)))+
    geom_jitter(data=cpo[!is.na(cpo.geno)&swarm=="A"], aes(x=temp.rack.cal, y=diapause.bin9, color=as.factor(cpo.geno)), height = 0.1, size=0.5)+
    labs(x="Temperature °C", y="", color="cpo genotype")+
    binomial_smooth()+
    scale_y_continuous(breaks=c(0,1))+
    scale_color_discrete(labels=c("A/A", "A/T", "T/T"))+
    theme(legend.position="none")


st10.B.cpo<-ggplot(data=cpo[!is.na(cpo.geno)&swarm=="B"], aes(x=temp.rack.cal, y=diapause.bin9, color=as.factor(cpo.geno)))+
    geom_jitter(data=cpo[!is.na(cpo.geno)&swarm=="B"], aes(x=temp.rack.cal, y=diapause.bin9, color=as.factor(cpo.geno)), height = 0.1, size=0.5)+
    labs(x="Temperature °C", y="", color="cpo geno")+
    binomial_smooth()+
    scale_y_continuous(breaks=c(0,1))+
    scale_color_discrete(labels=c("A/A", "A/T", "T/T"))+
    theme(legend.position=c(0.6, 0.7))



pdf("/scratch/pae3g/revisions/figures/cpo_tim.pdf", height=8, width=8)
plot_grid(st8.both.tim, st8.A.tim, st8.B.tim, st10.both.tim, st10.A.tim, st10.B.tim, st8.both.cpo, st8.A.cpo, st8.B.cpo, st10.both.cpo, st10.A.cpo, st10.B.cpo, nrow=4, labels=c("A", "", "", "B", "", "", "C", "", "","D", "", ""), rel_heights=c(.28, .24,.24, .24))
dev.off()


#stats for figure legend
summary(glm(diapause.bin9~final.geno+temp.rack.cal+photoperiod+swarm.y+generation+wolbachia, data=tim.genos, family=binomial))
summary(glm(diapause.bin9~final.geno+temp.rack.cal+generation+wolbachia, data=tim.genos[swarm.y=="A"], family=binomial))
summary(glm(diapause.bin9~final.geno+temp.rack.cal+generation+wolbachia, data=tim.genos[swarm.y=="B"], family=binomial)) #p=0.0189


summary(glm(diapause.bin~final.geno+temp.rack.cal+photoperiod+swarm.y+generation+wolbachia, data=tim.genos, family=binomial))
summary(glm(diapause.bin~final.geno+temp.rack.cal+generation+wolbachia, data=tim.genos[swarm.y=="A"], family=binomial))
summary(glm(diapause.bin~final.geno+temp.rack.cal+generation+wolbachia, data=tim.genos[swarm.y=="B"], family=binomial))


summary(glm(diapause.bin9~cpo.geno+temp.rack.cal+photoperiod+swarm+generation+wolbachia, data=cpo, family=binomial))
summary(glm(diapause.bin9~cpo.geno+temp.rack.cal+photoperiod+generation+wolbachia, data=cpo[swarm=="A"], family=binomial))
summary(glm(diapause.bin9~cpo.geno+temp.rack.cal+photoperiod+generation+wolbachia, data=cpo[swarm=="B"], family=binomial))


summary(glm(diapause.bin~cpo.geno+temp.rack.cal+photoperiod+swarm+generation+wolbachia, data=cpo, family=binomial))
summary(glm(diapause.bin~cpo.geno+temp.rack.cal+photoperiod+generation+wolbachia, data=cpo[swarm=="A"], family=binomial))
summary(glm(diapause.bin~cpo.geno+temp.rack.cal+photoperiod+generation+wolbachia, data=cpo[swarm=="B"], family=binomial))




#######################################
##### FIGURE S12: Machado 2019 ########
#######################################

load("/scratch/pae3g/revisions/evolution/bergland2019_sign_universal_threshold.Rdata") #on dryad
load("/scratch/pae3g/revisions/all_lasso_for_alan.Rdata") #on dryad

cline<-fread("/scratch/pae3g/revisions/evolution/east_coast_cline_V2_clean.txt") #on dryad
seas<-fread("/scratch/pae3g/revisions/evolution/seas_glm_NOTswitch_clean.txt")   #on dryad
freqs<-fread("/scratch/pae3g/oldscratch_recovered/evolution/east_coast_cline_V2_allele_freqs.txt") #on dryad

c<-merge(cline, seas, by=c("chr", "pos"))
c<-merge(c, freqs, by=c("chr", "pos"))

c<-c[f.hat>=0.05&f.hat<=0.95]
c<-c[clinal.beta<3&clinal.beta>(-3)]
c<-c[seas.beta<3&seas.beta>(-3)]





x<-merge(p, c[,.(chr, pos, clinal.beta, seas.beta, f.hat, clinal.p, seas.p)], by=c("chr", "pos"))

x[,ps.cline:=-1*coef*clinal.beta]
x[,ps.seas:=-1*coef*seas.beta]


x.sum<-x[,.(clinal=sum(ps.cline),
            seasonal=sum(ps.seas)),. (GRM, perm, pop, pheno)]

x.melt<-melt(x.sum, measure.vars=c("clinal", "seasonal"), id.vars=c("GRM", "perm", "pop", "pheno"))
setnames(x.melt, c("variable", "value", "GRM"), c("test", "poly" , "draw" ))
x.melt[,top:="lasso"]


b2019<-rbind(y[,.(poly, top, perm, draw, test, pheno, pop)], x.melt)
b2019[perm!=0&pop=="both", group:="both-permuted"]
b2019[perm==0&pop=="both", group:="both-observed"]
b2019[perm==0&pop=="A", group:="A-observed"]
b2019[perm==0&pop=="B", group:="B-observed"]
b2019[perm!=0&pop=="A", group:="A-permuted"]
b2019[perm!=0&pop=="B", group:="B-permuted"]

b2019[,pheno2:=ifelse(pheno=="diapause.bin", "st. 8", "st. 10")]
b2019[,pheno2:=factor(b2019$pheno2, levels=c("st. 8", "st. 10"))]


b2019.sum<-b2019[,.(med=median(poly), q.025=quantile(poly, 0.025), q.975=quantile(poly, .975)), .(group, pheno, pheno2, pop, test,top)]
b2019.sum[,pheno2:=factor(b2019.sum$pheno2, levels=c("st. 8", "st. 10"))]


b2019.sum[top==-5, th:="Top 0.001%"]
b2019.sum[top==-4, th:="Top 0.01%"]
b2019.sum[top==-3, th:="Top 0.1%"]
b2019.sum[top==-2, th:="Top 1%"]
b2019.sum[top==-1, th:="Top 10%"]
b2019.sum[top==0, th:="all SNPs"]
b2019.sum[top=="lasso", th:="LASSO"]
a.plot<-ggplot(b2019.sum[test=="clinal" & th%in%c("LASSO", "all SNPs", "Top 0.01%", "Top 0.1%")])+
    geom_point(data=b2019.sum[test=="clinal"&th%in%c("LASSO", "all SNPs", "Top 0.01%", "Top 0.1%")], aes(x=pop, y=med, color=group), position=position_dodge(width=0.5))+
    geom_errorbar(data=b2019.sum[test=="clinal"&th%in%c("LASSO", "all SNPs", "Top 0.01%", "Top 0.1%")], aes(x=pop, ymax=q.975, ymin=q.025, color=group), width=0.2, position=position_dodge(width=0.5))+
    labs(x="", y="sum(GWAS coefficient*model coefficient)", color="", title="clinal")+
    scale_color_manual(values=c("dodgerblue2", "grey80", "darkorchid", "grey80", "lightseagreen", "grey80"))+
    facet_grid(th~pheno2, scales ="free_y")+
    theme(legend.position = "none")

b.plot<-ggplot(b2019.sum[test=="seasonal" & th%in%c("LASSO", "all SNPs", "Top 0.01%", "Top 0.1%")])+
    geom_point(data=b2019.sum[test=="seasonal"&th%in%c("LASSO", "all SNPs", "Top 0.01%", "Top 0.1%")], aes(x=pop, y=med, color=group), position=position_dodge(width=0.5))+
    geom_errorbar(data=b2019.sum[test=="seasonal"&th%in%c("LASSO", "all SNPs", "Top 0.01%", "Top 0.1%")], aes(x=pop, ymax=q.975, ymin=q.025, color=group), width=0.2, position=position_dodge(width=0.5))+
    labs(x="", y="", color="", title="seasonal")+
    scale_color_manual(values=c("dodgerblue2", "grey80", "darkorchid", "grey80", "lightseagreen", "grey80"))+
    facet_grid(th~pheno2, scales ="free_y")+
    theme(legend.position = "none")

pdf("/scratch/pae3g/revisions/figures/machado2019_main.pdf", height=8, width=8)
plot_grid(a.plot, b.plot, labels=c("A", "B"), nrow=1)
dev.off()


##########################################################
##### FIGURE S13/S14 individual populations ##############
##########################################################



load("/scratch/pae3g/revisions/evolution/single_population_sign_universal_threshold.Rdata")

load("/scratch/pae3g/revisions/all_lasso_for_alan.Rdata")

load("/scratch/pae3g/oldscratch_recovered/evolution/core20delta.rdat")
pops<-names(deltas)

library(foreach)
pop.test<-foreach(pop=pops)%dopar%{
    print(pop)
    l<-merge(p, deltas[[pop]][is.finite(diff.logit)], by=c("chr", "pos"))
    #l[,TT:=ifelse(sign(coef)==1 & sign(diff.logit)==(1), T, F)]
    #l[,TF:=ifelse(sign(coef)==(-1) & sign(diff.logit)==1, T, F)]
    #l[,FT:=ifelse(sign(coef)==(1) & sign(diff.logit)==(-1), T, F)]
    #l[,FF:=ifelse(sign(coef)==(-1) & sign(diff.logit)==(-1), T, F)]
    l[,ps:=coef*diff.logit]
    return(l[,.(poly=sum(ps)), .(GRM, perm,population, pop, pheno)])
    
}
pop.test<-rbindlist(pop.test)
pop.test[,model:="nonloco"]
pop.test[,top:="lasso"]
setnames(pop.test, "GRM", "draw")
y[,or:=NULL]
x<-rbind(y, pop.test)

x[perm!=0&pop=="both", group:="both-permuted"]
x[perm==0&pop=="both", group:="both-observed"]
x[perm==0&pop=="A", group:="A-observed"]
x[perm==0&pop=="B", group:="B-observed"]
x[perm!=0&pop=="A", group:="A-permuted"]
x[perm!=0&pop=="B", group:="B-permuted"]



x.sum<-x[,.(med=median(poly, na.rm=T), q.025=quantile(poly, 0.025, na.rm=T), q.975=quantile(poly, .975, na.rm=T)), .( pheno, pop,top, population, group)]
x.sum[,pheno2:=ifelse(pheno=="diapause.bin","st. 8", "st. 10")]

x.sum[,pheno2:=factor(x.sum$pheno2, levels=c("st. 8", "st. 10"))]



x.sum[top==-5, th:="Top 0.001%"]
x.sum[top==-4, th:="Top 0.01%"]
x.sum[top==-3, th:="Top 0.1%"]
x.sum[top==-2, th:="Top 1%"]
x.sum[top==-1, th:="Top 10%"]
x.sum[top==0, th:="all SNPs"]
x.sum[top=="lasso", th:="LASSO"]

popinfo<-fread("/scratch/pae3g/revisions/popinfo_toplot.csv")
setnames(popinfo, "pop", "population")

x.sum<-merge(x.sum, popinfo, by="population")
x.sum[,label:=factor(x.sum$label, levels=popinfo[order(year)][order(latitude)][,label])]
a.plot<-ggplot(x.sum[pheno=="diapause.bin9" & th%in%c("LASSO", "all SNPs", "Top 0.01%", "Top 0.1%")])+
    geom_point(data=x.sum[pheno=="diapause.bin9"&th%in%c("LASSO", "all SNPs", "Top 0.01%", "Top 0.1%")], aes(x=pop, y=med, color=group), position=position_dodge(width=0.5))+
    geom_errorbar(data=x.sum[pheno=="diapause.bin9"&th%in%c("LASSO", "all SNPs", "Top 0.01%", "Top 0.1%")], aes(x=pop, ymax=q.975, ymin=q.025, color=group), width=0.2, position=position_dodge(width=0.5))+
    labs(x="", y="seasonal polygenic score", color="")+
    scale_color_manual(values=c("dodgerblue2", "grey80", "darkorchid", "grey80", "lightseagreen", "grey80"))+
    facet_grid(th~label, scales ="free_y")+
    theme(legend.position="none")

pdf("/scratch/pae3g/revisions/figures/individual_populations_stg10.pdf", height=8, width=20)
a.plot
dev.off()


b.plot<-ggplot(x.sum[pheno=="diapause.bin" & th%in%c("LASSO", "all SNPs", "Top 0.01%", "Top 0.1%")])+
    geom_point(data=x.sum[pheno=="diapause.bin"&th%in%c("LASSO", "all SNPs", "Top 0.01%", "Top 0.1%")], aes(x=pop, y=med, color=group), position=position_dodge(width=0.5))+
    geom_errorbar(data=x.sum[pheno=="diapause.bin"&th%in%c("LASSO", "all SNPs", "Top 0.01%", "Top 0.1%")], aes(x=pop, ymax=q.975, ymin=q.025, color=group), width=0.2, position=position_dodge(width=0.5))+
    labs(x="", y="seasonal polygenic score", color="")+
    scale_color_manual(values=c("dodgerblue2", "grey80", "darkorchid", "grey80", "lightseagreen", "grey80"))+
    facet_grid(th~label, scales ="free_y")+
    theme(legend.position="none")

pdf("/scratch/pae3g/revisions/figures/individual_populations_stg8.pdf", height=8, width=20)
b.plot
dev.off()




#################################
### FIGURE S16: admixture #######
#################################


load("/scratch/pae3g/revisions/evolution/ZI_admix_universal_threshold.Rdata")

admix.th<-copy(admix)
load("/scratch/pae3g/revisions/evolution/admix_universal_lasso.Rdata")
admix[,top:="lasso"]
admix.sum<-admix[,.(total.admix=sum(n.admix), unique.admix=sum(n.admix>0), ntotal=.N), .(perm, draw, pheno, pop, top)]

admix<-rbind(admix.sum, admix.th[,.(perm, draw, pheno, pop, top, unique.admix, total.admix, ntotal)])

admix[,ratio:=total.admix/ntotal]
admix[,prop:=unique.admix/ntotal]

admix[perm!=0&pop=="both", group:="both-permuted"]
admix[perm==0&pop=="both", group:="both-observed"]
admix[perm==0&pop=="A", group:="A-observed"]
admix[perm==0&pop=="B", group:="B-observed"]
admix[perm!=0&pop=="A", group:="A-permuted"]
admix[perm!=0&pop=="B", group:="B-permuted"]



admix.melt<-melt(admix, measure.vars=c('ratio', 'prop'), id.vars=c("perm", "draw", "pheno", "pop", 'top', 'group'))
admix.melt.sum<-admix.melt[,.(med=median(value, na.rm=T), q.025=quantile(value, 0.025, na.rm=T), q.975=quantile(value, .975, na.rm=T)), .( pheno, pop,top, variable, group)]
admix.melt.sum[,pheno2:=ifelse(pheno=="diapause.bin","st. 8", "st. 10")]

admix.melt.sum[,pheno2:=factor(admix.melt.sum$pheno2, levels=c("st. 8", "st. 10"))]



admix.melt.sum[top==-5, th:="Top 0.001%"]
admix.melt.sum[top==-4, th:="Top 0.01%"]
admix.melt.sum[top==-3, th:="Top 0.1%"]
admix.melt.sum[top==-2, th:="Top 1%"]
admix.melt.sum[top==-1, th:="Top 10%"]
admix.melt.sum[top==0, th:="all SNPs"]
admix.melt.sum[top=="lasso", th:="LASSO"]


a.plot<-ggplot(admix.melt.sum[variable=="prop" & th%in%c("LASSO", "all SNPs", "Top 0.01%", "Top 0.1%")])+
    geom_point(data=admix.melt.sum[variable=="prop"&th%in%c("LASSO", "all SNPs", "Top 0.01%", "Top 0.1%")], aes(x=pop, y=med, color=group), position=position_dodge(width=0.5))+
    geom_errorbar(data=admix.melt.sum[variable=="prop"&th%in%c("LASSO", "all SNPs", "Top 0.01%", "Top 0.1%")], aes(x=pop, ymax=q.975, ymin=q.025, color=group), width=0.2, position=position_dodge(width=0.5))+
    labs(x="", y="proportion of SNPs in admixture tracts", color="")+
    scale_color_manual(values=c("dodgerblue2", "grey80", "darkorchid", "grey80", "lightseagreen", "grey80"))+
    facet_grid(th~pheno2, scales ="free_y")+
    theme(legend.position="none")

pdf("/scratch/pae3g/revisions/figures/admixture_proportion.pdf", height=6, width=5)
a.plot
dev.off()



################################################
############# FIGURE S18 #######################
################################################


keep.paths<-fread( "/scratch/pae3g/evolution/all_reconstructed_haplotypes_melted.txt")
prop.line<-keep.paths[, .(line.length=sum(path.length)),.(line, chromosome, swarm, generation)]
prop.line[,proportion:=line.length/sum(as.numeric(line.length[chromosome==chromosome]))]

#cacluate total length of each chromosome in each swarm and add to file
for (chr in c("2L", "2R", "3L", "3R", "X")) {
    for (cage in c("A", "B")){
        for(gen in c(4,5)){
            chr.total=sum(prop.line[chromosome==chr&swarm==cage&generation==gen, line.length], na.rm=TRUE)
            prop.line[chromosome==chr&swarm==cage&generation==gen, percent:=line.length/chr.total]
        }
    }
}


#read in geography and reorder it as a factor

geo<-fread("/scratch/pae3g/evolution/strain_geography.csv") #available in Table S2
setnames(geo, "strain", "line")
geo[,geography:=factor(geography, levels=c("Carribean", "Southeast", "DGRP", "PA-Fall", "PA-Spring", "Ithica", "Maine"))]
prop.line<-merge(prop.line, geo, by="line")

#order line names by geography
sorted_lines<-geo[order(geography), line]
prop.line[,line:=factor(line, levels=sorted_lines)]
prop.line[,geography:=factor(geography, levels=c("Carribean", "Southeast", "DGRP", "PA-Fall", "PA-Spring", "Ithica", "Maine"))]

#make factor for cage+ generation
prop.line[,group.name:=paste("cage ",swarm, " gen. ", generation)]


#plot with F4s and F5s mirrored above and below the axis

pdf("/scratch/pae3g/prop_each_line.pdf", height=8, width=8)
ggplot(prop.line[generation==4], aes(x=line, y=percent*100, fill=geography))+
    geom_bar(stat="identity")+
    geom_bar(data=prop.line[generation==5], aes(x=line, y=percent*-100, fill=geography), stat="identity")+
    facet_grid(chromosome~swarm, scales="free_x")+
    theme(axis.text.x=element_text(angle=45,hjust=1, size=8))+
    labs(y="% of chromosome derived from founder line", fill="Founder origin")+
    geom_hline(yintercept=1/34*100, color="grey", linetype="dashed")+
    geom_hline(yintercept=0, color="black")+
    geom_hline(yintercept=-1/34*100, color="grey", linetype="dashed")+
    scale_fill_discrete(labels=c("Carribean", "Southeast", "North Carolina", "Penn-Fall", "Penn-Spring", "New York", "Maine"))

dev.off()


################################################
########### FIGURE S19 #########################
################################################

#this file needs to be collapsed to create consecutive paths
a <- fread("/scratch/pae3g/genome-reconstruction/simulated_actual_haplotypes.txt")

a[,run:=rleid(id, chromosome, haplotype, lineID, swarm, rep, gen)]
all.paths.cons=a[,.(cons.start=min(start), cons.stop=max(stop)), .(id, chromosome, haplotype, lineID, swarm, run, rep, gen)]
all.paths.cons[,path.length:=cons.stop-cons.start]

sim<-all.paths.cons

#melt these into consecutive paths
all.sim.rec.melts<-fread("/scratch/pae3g/genome-reconstruction/all_simulated_reconstructed_haplotypes_melted.txt")
all.sim.rec.cons=all.sim.rec.melt[,.(cons.start=min(start), cons.stop=max(stop)), .(sample.id, chromosome, haplotype, line, swarm, run, gen)]
all.sim.rec.cons[,path.length:=cons.stop-cons.start]


#clean this up same was as actual data

all.sim.rec.cons[path.length<1000000, line:="UNKNOWN"]
all.sim.rec.cons<-all.sim.rec.cons[order(sample.id, chromosome, haplotype, cons.start)]
all.sim.rec.cons[,run:=rleid(sample.id, chromosome, haplotype, line, gen)]

all.sim.rec.cons.merge=all.sim.rec.cons[,.(cons.start=min(cons.start), cons.stop=max(cons.stop)), .(sample.id, chromosome, haplotype, line, swarm, run,  gen)]

#now "UNKNOWNS" are bridged into single paths
all.sim.rec.cons.merge[,path.length:=cons.stop-cons.start]

#drop singleton short paths and merge again

all.sim.rec.cons.merge<-all.sim.rec.cons.merge[path.length>1000000]
all.sim.rec.cons.merge[,run:=rleid(sample.id, chromosome, haplotype, line)]

all.sim.rec.cons.merge2=all.sim.rec.cons.merge[,.(cons.start=min(cons.start), cons.stop=max(cons.stop)), .(sample.id, chromosome, haplotype, line, swarm, run,  gen)]

all.sim.rec.cons.merge2[,path.length:=cons.stop-cons.start]
#read in actual paths

real<-fread("/scratch/pae3g/evolution/all_reconstructed_haplotypes_melted.txt")
setnames(real, "generation", "gen")

#clean up "real" data by merging and bridging over short paths
real[path.length<1000000, line:="UNKNOWN"]
real<-real[order(sample.id, chromosome, haplotype, cons.start)]
real[,run:=rleid(sample.id, chromosome, haplotype, line)]

real.merge=real[,.(cons.start=min(cons.start), cons.stop=max(cons.stop)), .(sample.id, chromosome, haplotype, line, swarm, run,  gen)]

#now "UNKNOWNS" are bridged into single paths
real.merge[,path.length:=cons.stop-cons.start]

#drop singleton short paths and merge again

real.merge<-real.merge[path.length>1000000]
real.merge[,run:=rleid(sample.id, chromosome, haplotype, line)]

real.merge2=real.merge[,.(cons.start=min(cons.start), cons.stop=max(cons.stop)), .(sample.id, chromosome, haplotype, line, swarm, run,  gen)]

real.merge2[,path.length:=cons.stop-cons.start]

#count recombinations in each dataset
sim.r<-sim[,.(R=.N-1),.(id,chromosome, swarm, gen)]
sim.rec.r<-all.sim.rec.melt[,.(R=.N-1), .(sample.id, chromosome, swarm, gen)]
sim.rec.merge2.r<-all.sim.rec.cons.merge2[,.(R=.N-1), .(sample.id, chromosome, swarm, gen)]
real.r<-real[,.(R=.N-1),.(sample.id,chromosome, swarm, gen)]
real.merge2.r<-real.merge2[,.(R=.N-1),.(sample.id,chromosome, swarm, gen)]


#real.merge2 has short paths merged into longer unknowns and "bridges" drawn where a short path falls in the middle of two paths of the same parent

plot1<-ggplot(sim, aes(x=log10(path.length), color="Simulation"))+
    stat_density(size=1.5, geom="line")+facet_wrap(~gen)+
    stat_density(data=real, aes(x=log10(path.length), color="Reconstruction"), geom="line")+
    stat_density(data=real.merge2, aes(x=log10(path.length), color="Cleaned up"), geom="line")+
    labs(x=expression("log"[10]*"(Haplotype length, bp)"), y="density", linetype="Population", color=NULL, title="Actual hybrid swarms")+
    scale_color_manual(values=c("#F8766D", "#00BFC4", "grey"))+
    theme(legend.position = "none")+
    lims(x=c(4,8))

plot2<-ggplot(sim, aes(x=log10(path.length), color="Simulation"))+
    stat_density(size=1.5,  geom="line")+facet_wrap(~gen)+
    stat_density(data=all.sim.rec.melt, aes(x=log10(path.length), color="Reconstruction"), geom="line")+
    stat_density(data=all.sim.rec.cons.merge2, aes(x=log10(path.length), color="Cleaned up"), geom="line")+
    labs(x=expression("log"[10]*"(Haplotype length, bp)"), y="density", linetype="Population", color=NULL, title="Simulated hybrid swarms")+
    scale_color_manual(values=c("#7CAE00", "#C77CFF", "grey"))+
    theme(legend.position = "none")+
    lims(x=c(4,8))


plot3<-ggplot(sim.r, aes(x=R, color="Simulation"))+
    stat_density(size=1.5, adjust=6, geom="line")+facet_wrap(~gen)+
    stat_density(data=real.r, aes(x=R, color="Reconstruction"), geom="line", adjust=6)+
    stat_density(data=real.merge2.r, aes(x=R, color="Cleaned up"), adjust=6 ,geom="line")+
    labs(x="Recombinations per diploid chromosome", y="density", linetype="Population", color=NULL)+
    scale_color_manual(values=c("#F8766D", "#00BFC4", "grey"))+
    theme(legend.position = "bottom", legend.text=element_text(size=9))+
    lims(x=c(0,20))

plot4<-ggplot(sim.r, aes(x=R, color="Simulation"))+
    stat_density(size=1.5, adjust=6,  geom="line")+facet_wrap(~gen)+
    stat_density(data=sim.rec.r, aes(x=R, color="Reconstruction"), adjust=6, geom="line")+
    stat_density(data=sim.rec.merge2.r, aes(x=R, color="Cleaned up"), adjust=6, geom="line")+
    scale_color_manual(values=c("#7CAE00", "#C77CFF", "grey"))+
    labs(x="Recombinations per diploid chromosome", y="density", linetype="Population", color=NULL)+
    theme(legend.position = "bottom", legend.text=element_text(size=9))+
    lims(x=c(0,20))


l<-plot_grid(plot1, plot3, nrow=2, labels=c("A", "C"))
r<-plot_grid(plot2, plot4, nrow=2, labels=c("B", "D"))
pdf("/scratch/pae3g/FigureS5.pdf", height=7.5, width=8.75)

plot_grid(l,r,nrow=1)
dev.off()

#K-S test on distributions

#haplotype length real original vs sim
ks.test(real[gen==4,path.length], sim[gen==4, path.length]) #P< 2x10-16, D=0.2
ks.test(real[gen==5,path.length], sim[gen==5, path.length])#P< 2x10-16 D=0.21

#haplotype length real cleaned vs sim
ks.test(real.merge2[gen==4,path.length], sim[gen==4, path.length]) #P< 2x10-16, D=.07
ks.test(real.merge2[gen==5,path.length], sim[gen==5, path.length])#P< 2x10-16, D=.05

# recomb real original vs sim
ks.test(real.r[gen==4,R], sim.r[gen==4, R]) #P< 2x10-16, D=0.1
ks.test(real.r[gen==5,R], sim.r[gen==5, R])#P< 2x10-16 D=0.13

# reeomb real cleaned vs sim
ks.test(real.merge2.r[gen==4,R], sim.r[gen==4, R]) #P< 6 x 10-10, D=.06
ks.test(real.merge2.r[gen==5,R], sim.r[gen==5, R])#P< 8 x 10-6, D=.05


#haplotype length simulated original vs sim
ks.test(all.sim.rec.melt[gen==4,path.length], sim[gen==4, path.length]) #P< 2x10-16, D=0.2
ks.test(all.sim.rec.melt[gen==5,path.length], sim[gen==5, path.length])#P< 2x10-16 D=0.19

#haplotype length simulated cleaned vs sim
ks.test(all.sim.rec.cons.merge2[gen==4,path.length], sim[gen==4, path.length]) #P< 2x10-16, D=.065
ks.test(all.sim.rec.cons.merge2[gen==5,path.length], sim[gen==5, path.length])#P< 2x10-16, D=.05

# recomb real simulated vs sim
ks.test(sim.rec.r[gen==4,R], sim.r[gen==4, R]) #P< 2x10-16, D=0.34
ks.test(sim.rec.r[gen==5,R], sim.r[gen==5, R])#P< 2x10-16 D=0.40

# recomb real simulated vs sim
ks.test(sim.rec.merge2.r[gen==4,R], sim.r[gen==4, R]) #P=.09, D=.03
ks.test(sim.rec.merge2.r[gen==5,R], sim.r[gen==5, R])#P< .33 D=.02

#################################  
### FIGURE S20: accuracy ########  
################################  


#on rivanna

library(data.table)
library(foreach)
library(cowplot)
library(ggbeeswarm)

acc<-fread("/scratch/pae3g/genome-reconstruction/accuracy.txt")


acc.sum<-acc[, .(sites=sum(nSites), matches=sum(nMatches.14)), .(ind_id, group)]
acc.sum[,accuracy:=matches/sites]

acc[,.(med=median(percentMatches.14), q=quantile(percentMatches.14, .05)), .(group)]

pdf("/scratch/pae3g/FigureS6.pdf", height=6, width=6)
ggplot(acc.sum, aes(x=group, y=accuracy))+
    geom_beeswarm(cex= 0.5, alpha=0.5)+
    scale_y_continuous(trans='logit', breaks=c(0.5, 0.75, 0.9, 0.99, 0.999, 0.9999, 0.99999))+
    labs(y="Proportion of correctly reconstructed genotypes", 
         x="Simulated Population and Generation")

dev.off()






##############################
##### FIGURE S21: IBS ########
##############################

library(GWASTools)
library(data.table)
library(SNPRelate)
library(gdsfmt)
library(foreach)
library(viridis)
library(cowplot)
geno<-snpgdsOpen("/scratch/pae3g/genome-reconstruction/final3.vcf.gds")


#read genotypes and filters
filters=fread("/scratch/pae3g/final_reconstruction2/hwe_missing_maf_filters.txt")

maf <- 0.05
missing.rate <- 0.15
threads <- 20
pass=filters[qc_filter=="PASS", snp.id]

phenos <- fread("/scratch/pae3g/phenos/phenos_062018.txt")
setkey(phenos, sample.id)

genotyped_samples<- read.gdsn(index.gdsn(geno, "sample.id"))
#only use genotypes from swarm

ids.to.use=intersect(genotyped_samples, phenos$sample.id)

snpset <- snpgdsLDpruning(geno, 
                          ld.threshold=0.2, 
                          slide.max.bp = 5000, 
                          autosome.only=FALSE,
                          missing.rate=.15,
                          sample.id=ids.to.use,
                          snp.id=pass,
                          maf=.05)


snpset.use <- unlist(snpset)

a<- snpgdsIBS(geno, snp.id=snpset.use)
b<-a$ibs
rownames(b)=a$sample.id
colnames(b)=a$sample.id

#phenos<-merge(phenos, pca.dt, by="id")
#sample.order<-phenos[id%in%genotyped_samples][order(pc1)][order(swarm)][,id]

a.melt<-data.table(melt(b))
setnames(a.melt, "Var1", "id")
a.melt<-merge(a.melt, phenos[,.(id, swarm)], by="id")

sample.order<-unique(a.melt[order(value)][order(swarm)]$id)

a.melt[,id:=factor(a.melt$id, levels=sample.order)]
a.melt[,Var2:=factor(a.melt$Var2, levels=sample.order)]

tiff("/scratch/pae3g/IBS_GRM.tif", height=6, width=6, units="in", res=300, compression = "lzw")
ggplot(a.melt, aes(x=id, y=Var2, fill=value))+
    geom_tile()+
    scale_fill_viridis(option="viridis")+ 
    #scale_fill_gradientn(colours = terrain.colors(20))+
    theme(line = element_blank(),
          axis.text = element_blank(),
          title = element_blank())+
    labs(fill="IBS")

dev.off()





