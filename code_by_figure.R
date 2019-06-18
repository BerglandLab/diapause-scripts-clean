#Erickson et al 2019 figures


#################################
#### FIGURE 1 ###################
#################################

#make a map of starting lines
library(maps)
library(data.table)
library(cowplot)

#these data are available in File S2

pops<-fread("~/Box Sync/HSparents/strain_geography.csv")
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

p<-fread("~/Box Sync/hybridSwarm/phenos_062018.txt")
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

#pull in my diagram of ovary development
a <- ggdraw() + draw_image("~/Box Sync/manuscripts/egg_stage_color_coded.jpg", scale = 0.9)

#plot ovary phenotype as a function of temperature group
b<-ggplot(p[!is.na(diapause.group)], aes(temp.group))+geom_bar(aes(fill=diapause.group))+
    labs(x="Temperature", 
         fill="Ovary\nPhenotype", 
         y="Number of Flies", 
         title="All flies")+
    theme(axis.text.x=element_text(angle=45,hjust=1))

#plot ovary stages of individuals with eggs
c<-ggplot(p[!is.na(egg.group)], aes(temp.group))+geom_bar(aes(fill=egg.group))+
    labs(x="Temperature", 
         fill="Ovary\nPhenotype", 
         y="Number of Flies", 
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

#plot models of likelihood of diapause as a function of temperature for different photoperiods
d<-ggplot(p, aes(x=temp.rack.cal, y=diapause.bin, color=as.factor(photoperiod)))+
    binomial_smooth()+
    labs(x="Temperature °C", 
         y="Probability of diapause", 
         color="Photoperiod",
         title="Stage 8")+
    lims(y=c(0,1))+
    geom_point(data=p, aes(x=temp.rack.cal, y=diapause.bin), color="grey", alpha=0.1)


e<-ggplot(p, aes(x=temp.rack.cal, y=diapause.bin9, color=as.factor(photoperiod)))+
    binomial_smooth()+
    labs(x="Temperature °C", 
         y="Probability of diapause", 
         color="Photoperiod",
         title="Stage 10")+
    lims(y=c(0,1))+
    geom_point(data=p, aes(x=temp.rack.cal, y=diapause.bin), color="grey", alpha=0.1)
f<-ggplot(p, aes(x=temp.rack.cal, y=eggP, color=as.factor(photoperiod)))+
    binomial_smooth()+
    labs(x="Temperature °C",
         y="Probability of diapause", 
         color="Photoperiod",
         title="Stage 14")+
    lims(y=c(0,1))+
    geom_point(data=p, aes(x=temp.rack.cal, y=diapause.bin), color="grey", alpha=0.1)

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
pdf("~/Desktop/figure2.pdf", height=6, width=8)
plot_grid(top.row, bottom.row, nrow=2)
dev.off()


#stats for figure legend for this and supplemental figures

summary(aov(glm(diapause.bin~temp.rack.cal+(photoperiod)+swarm+generation, data=p, family="binomial")))
aov(glm(diapause.bin~temp.rack.cal+(photoperiod)+swarm+generation, data=p, family="binomial"))
summary(aov(glm(diapause.bin9~temp.rack.cal+(photoperiod)+swarm+generation, data=p, family="binomial")))
summary(aov(glm(eggP~temp.rack.cal+(photoperiod)+swarm+generation, data=p, family="binomial")))


############################
######## FIGURE 3 #########
############################




## 16 core largemem for lambda/qq stuff
library(data.table)
library(foreach)
library(doMC)
registerDoMC(16)
library(cowplot)
library(ggbeeswarm)

gwas<-fread("/scratch/pae3g/evolution/gwas_p_score_inv_id_perm0.txt.txt")

#adjust p values and add groupings for color coding the plt
gwas[,color:=1]
gwas[chr=="2R"|chr=="3R", color:=2]
gwas[,p.adj:=p.adjust(gwas.p, method="fdr")]
gwas[p.adj<.05, fdr.group:="FDR < .05"]
gwas[p.adj<.025, fdr.group:="FDR < .025"]
gwas[,q:=frank(gwas.p)/length(gwas$gwas.p)]
gwas[p.adj>.05&q<0.005, fdr.group:="top 0.5%"]

#determine locations to draw breaks by finding the midpoint of each chromosome
chr.table<-gwas[,.(max.id=max(id)), .(chr)]
chr.table<-foreach(i=1:5, .combine="rbind")%do%{
    temp<-chr.table[i]
    temp[,midpoint:=(chr.table[i,max.id]+chr.table[i-1, max.id])/2]
    return(temp)
}
chr.table[chr=="2L", midpoint:=max.id/2]
chr.table[,fdr.group:=NA]


# FINAL PLOT


mh<-ggplot(gwas[is.na(fdr.group)|fdr.group=="top 0.5%"], aes(x=id, y=-log10(gwas.p), color=as.factor(color)))+
    geom_point()+
    scale_color_manual(values=c("black", "grey40", "#7CAE00", "#00BFC4", "#C77CFF", "#F8766D"))+
    geom_point(data=gwas[fdr.group=="FDR < .025" | fdr.group=="FDR < .05"], aes(x=id, y=-log10(gwas.p), color=fdr.group))+
    theme(legend.position="none")+labs(x=NULL, y=expression("-log"[10]*"(P)"))+scale_x_continuous(breaks=chr.table$midpoint, labels=c("2L", "2R", "3L", "3R", "X"))

##qqplot of permutation

dat<-foreach(perm=c(0,101:200))%dopar%{
    print(perm)
    d<-fread(paste0("/scratch/pae3g/evolution/gwas_p_score_inv_id_perm", perm, ".txt"))
    d[,perm:=perm]
    d[,o := -log10(sort(gwas.p,decreasing=F))]
    d[,e := -log10( 1:length(gwas.p)/length(gwas.p))]
    return(d)
}

dat<-rbindlist(dat)

#subsample 1/10000th of data with p>0.05 for qqplots

dat.sig<-dat[perm!=0&o>(-log10(0.05))]
dat.low<-dat[perm!=0&o<(-log10(0.05))]

dat.low<-dat.low[sample(c(1:nrow(dat.low)),nrow(dat.low)/1000)]

dat.perm<-rbind(dat.sig, dat.low)
#replace infinities

qq<-ggplot(dat[perm==0], aes(x=e, y=o) )+
    geom_point(color="lightseagreen")+
    labs(x=expression("-log"[10]*"(expected P)"), 
         y=expression("-log"[10]*"(observed P)"))+
    geom_point(data=dat.perm, aes(x=e, y=o),  color="grey", size=0.1)+
    geom_abline(slope=1, intercept=0)



##calculate lambda gc values for gwas permutations


y<-fread("/scratch/pae3g/evolution/lambda.txt")
h<-fread("/scratch/pae3g/evolution/gcta_hsq.txt")


y[,Source:="λ GC"]
lambda.plot<-ggplot(y[perm!=0], aes(x=Source, y=lambda))+
    geom_quasirandom(color="grey40", varwidth = TRUE, size=.5, method = "smiley")+
    geom_point(data=y[perm==0], aes(x=Source, y=lambda), color="lightseagreen", size=2)+labs(y=NULL, x=NULL)

h.plot<-ggplot(h[perm!=0&(Source=="V(G)"|Source=="V(e)"|Source=="V(G)/Vp")], aes(x=Source, y=Variance))+
    geom_quasirandom(color="grey40", varwidth = TRUE, size=.5, method = "smiley")+
    geom_point(data=h[perm==0&(Source=="V(G)"|Source=="V(e)"|Source=="V(G)/Vp")], aes(x=Source, y=Variance), color="lightseagreen", size=2)+labs(y=NULL, x=NULL)+scale_x_discrete(labels=c(expression("V"["e"]), expression("V"["g"]), expression("h"^2)))

bottom<-plot_grid(qq, lambda.plot, h.plot, rel_widths =c(0.5, 0.15, 0.35), labels = c("B", "C", "D"), nrow=1, align="h")



jpeg("/scratch/pae3g/Figure3.jpeg", height=8, width=7.75, res = 300, units = "in")
plot_grid(mh, bottom, nrow=2, labels=c("A",""), rel_heights = c(.5, .5))
dev.off()


####################################
##### FIGURE 4 #####################
####################################

#LD heatmaps of top SNPs

library(data.table)
library(SNPRelate)
library(foreach)
library(GWASTools)
library(cowplot)
library(viridis)

#open an example imputed genotype file
geno <- snpgdsOpen("/scratch/pae3g/genome-reconstruction/final2_draw1_replaced.vcf.gds", allow.fork=T)

#pull out snp info from this file
a<-snpgdsSNPList(geno)
info<-data.table(chr=a$chromosome,
                 pos=a$position,
                 snp.id=a$snp.id)
gwas<-fread("/scratch/pae3g/evolution/gwas_p_score_inv_id.txt")

gwas<-merge(gwas, info, by=c("chr", "pos"))
gwas[,fdr:=p.adjust(gwas.p, method="fdr")]


top.025<-snpgdsLDMat(geno, snp.id=gwas[fdr<.025, snp.id], slide=-1)

a.melt<-data.table(melt(top.025$LD))
a<-ggplot(a.melt, aes(x=Var1, y=Var2, fill=value*value))+
    geom_tile()+
    scale_fill_viridis(option="viridis")+ 
    theme(line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank())+
    labs(title="FDR < 0.025", x="", y="", fill=expression("R"^2))


top.05<-snpgdsLDMat(geno, snp.id=gwas[fdr<.05, snp.id], slide=-1)

b.melt<-data.table(melt(top.05$LD))
b<-ggplot(b.melt, aes(x=Var1, y=Var2, fill=value*value))+
    geom_tile()+
    scale_fill_viridis(option="viridis")+ 
    theme(line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank())+
    labs(title="FDR < 0.05", x="", y="")


leg<-get_legend(a)

ld_heatmaps<-plot_grid(a+theme(legend.position="none"),
                       b+theme(legend.position="none"),
                       leg,
                       nrow=1,
                       rel_widths = c(.4,.4,.1),
                       labels=c("E", "F"))

#plot minor allele frquency spectra for DGRP and hybrid swarm

#using genotypes open above

swarm.freqs<-snpgdsSNPRateFreq(geno, snp.id=gwas$snp.id)

snpgdsClose(geno)

geno<-snpgdsOpen("/scratch/pae3g/evolution/dgrp2.vcf.gds", allow.fork=T)
dgrp.freqs<-snpgdsSNPRateFreq(geno, snp.id=gwas$snp.id)

dat<-data.table(freq=c(swarm.freqs$MinorFreq, dgrp.freqs$MinorFreq),
                id=c(rep("hybrid swarm", length(swarm.freqs$MinorFreq)), 
                     rep("DGRP", length(dgrp.freqs$MinorFreq))))

maf.hist<-ggplot(data=dat[id=="DGRP"], aes(x=freq, fill=id))+geom_histogram(binwidth=0.01, alpha=0.5, aes(y = stat(count / sum(count))))+labs(fill=NULL, x="Minor Allele Frequency", y="Proportion", title=" ")+geom_histogram(data=dat[id=="hybrid swarm"], aes(x=freq, fill=id, y = stat(count / sum(count))), binwidth=0.01, alpha=0.5)+ theme(legend.position = c(0.3, 0.5))


row2<-plot_grid(maf.hist, ld_heatmaps, rel_widths = c(0.3, 0.7), labels=c("D",""))

#index SNPs 

library(data.table)
library(foreach)
library(cowplot)

blocks<-foreach(ld=c("0.2", "0.3", "0.4", "0.5", "0.6"))%do%{
    b<-foreach(pval=c("0.000001065899", "0.00001094209", "0.0002394024"))%do%{
        c<-foreach(kb=c(10, 25, 50,100,200,300,400,500))%do%{
            blocks=fread(paste0("/scratch/pae3g/haplotypes/clumped_p", pval, "_ld", ld, "_kb", kb, ".clumped"))
            blocks[,pval:=pval]
            blocks[,ld:=ld]
            blocks[,kb:=kb]
            return(blocks)
        }
        return(rbindlist(c))
    }
    return(rbindlist(b))
}

blocks<-rbindlist(blocks)
blocks[pval=="0.000001065899", fdr:=0.025]
blocks[pval=="0.00001094209", fdr:=0.05]
blocks[pval=="0.0002394024", fdr:=0.1]


blocks.sum<-blocks[,.(n.clumps=.N), .(ld, pval, kb, fdr)]


e<-ggplot(blocks.sum[fdr==.025], aes(x=kb, y=n.clumps, color=ld, group=ld))+geom_point()+geom_line()+labs(x="max. distance, kb", y="#index SNPs", title="FDR < 0.025, n = 101", color="LD\nthreshold")

f<-ggplot(blocks.sum[fdr==.05], aes(x=kb, y=n.clumps, color=ld, group=ld))+geom_point()+geom_line()+labs(x="max. distance, kb", y="", title="FDR < 0.05, n = 490", color="LD\nthreshold")

g<-ggplot(blocks.sum[fdr==.1], aes(x=kb, y=n.clumps, color=ld, group=ld))+geom_point()+geom_line()+labs(x="max. distance, kb", y="", title="FDR < 0.1, n = 5316", color="LD\nthreshold")


l<-get_legend(e)


row3<-plot_grid(e+theme(legend.position="none"),
                f+theme(legend.position="none"),
                g+theme(legend.position="none"),
                l,
                nrow=1,
                rel_widths = c(0.3,0.3,0.3, 0.1),
                labels=c("G", "H", "I"))


#plot ld decay

s.c<-fread("/scratch/pae3g/evolution/ld_decay_maf0.1_byswarm.txt")
d.c<-fread("/scratch/pae3g/evolution/ld_decay_dgrp_maf0.1.txt")
d.c[,rsq:=ld*ld]
s.c[,rsq:=ld*ld]

s.c.sum<-s.c[,.(med.rsq=median(rsq, na.rm=T)), .(dist)]
d.c.sum<-d.c[,.(med.rsq=median(rsq, na.rm=T)), .(dist)]

common<-ggplot(d.c.sum[dist>50], aes(x=log10(dist), y=log10(med.rsq), color="DGRP"))+
    geom_point()+
    geom_line()+
    geom_point(data=s.c.sum[dist>50], aes(x=log10(dist), y=log10(med.rsq), color="hybrid swarm"))+
    geom_line(data=s.c.sum[dist>50], aes(x=log10(dist), y=log10(med.rsq), color="hybrid swarm"))+
    labs(x=expression("log"[10]*"(distance, bp)"), y=" ", title="MAF > 0.1", color=NULL)+
    lims(y=c(-4, -.5))


s.r<-fread("/scratch/pae3g/evolution/ld_decay_mafunder.05_byswarm.txt")
d.r<-fread("/scratch/pae3g/evolution/ld_decay_dgrp_mafunder.05.txt")
d.r[,rsq:=ld*ld]
s.r[,rsq:=ld*ld]

s.r.sum<-s.r[,.(med.rsq=median(rsq, na.rm=T)), .(dist)]
d.r.sum<-d.r[,.(med.rsq=median(rsq, na.rm=T)), .(dist)]

rare<-ggplot(d.r.sum[dist>50], aes(x=log10(dist), y=log10(med.rsq), color="DGRP"))+
    geom_point()+
    geom_line()+
    geom_point(data=s.r.sum[dist>50], aes(x=log10(dist), y=log10(med.rsq), color="hybrid\nswarm"))+
    geom_line(data=s.r.sum[dist>50], aes(x=log10(dist), y=log10(med.rsq), color="hybrid\nswarm"))+
    labs(x=expression("log"[10]*"(distance, bp)"), y=" ", title="MAF < 0.05", color=NULL)+
    lims(y=c(-4, -.5))

s.a<-fread("/scratch/pae3g/evolution/ld_decay_all_byswarm.txt")
d.a<-fread("/scratch/pae3g/evolution/ld_decay_dgrp_all.txt")
d.a[,rsq:=ld*ld]
s.a[,rsq:=ld*ld]

s.a.sum<-s.a[,.(med.rsq=median(rsq, na.rm=T)), .(dist)]
d.a.sum<-d.a[,.(med.rsq=median(rsq, na.rm=T)), .(dist)]

all<-ggplot(d.a.sum[dist>50], aes(x=log10(dist), y=log10(med.rsq), color="DGRP"))+
    geom_point()+
    geom_line()+
    geom_point(data=s.a.sum[dist>50], aes(x=log10(dist), y=log10(med.rsq), color="hybrid swarm"))+
    geom_line(data=s.a.sum[dist>50], aes(x=log10(dist), y=log10(med.rsq), color="hybrid swarm"))+
    labs(x=expression("log"[10]*"(distance, bp)"), y=expression("log"[10]*"(median R"^2*")"), title="All SNPs", color=NULL)+
    lims(y=c(-4, -.5))

le<-get_legend(rare)
row1<-plot_grid(all+theme(legend.position="none"),
                rare+theme(legend.position="none"),
                common+theme(legend.position="none"),
                le,
                nrow=1,
                labels=c("A", "B", "C", ""),
                rel_widths = c(0.3, 0.3, 0.3, 0.1))
pdf("/scratch/pae3g/Figure6.pdf", height=8, width=8)
plot_grid(row1,row2, row3, nrow=3)
dev.off()



###################
### FIGURE 5 ######
###################
library(data.table)
library(cowplot)
library(viridis)
library(ggbeeswarm)


clump.enrich<-fread("/scratch/pae3g/evolution/bergland2014_enrich_intersectfirst_clump.txt")

#since the number of index SNPs is different in each permutation, we need to calculate the fraction of sites that are clinal/seasonal rather than just the number

clump.enrich[,prop.test:=TT/(TT+TF)]

#summarize all permutations except 0, which is original ordering of data
perm.sum<-clump.enrich[perm!=0, .(med.TT=median(TT, na.rm=T),
                                  med.prop.test=median(prop.test, na.rm=T),
                                  med.prop.conc=median(conc.actual, na.rm=T),
                                  q.95.TT=quantile(TT, .95, na.rm=T),
                                  q.95.prop=quantile(prop.test, .95, na.rm=T),
                                  q.95.conc=quantile(conc.actual, .95, na.rm=T)),
                       .(clinal.th, gwas.th, seas.th, test)]

#merge with actual data
perm.sum<-merge(perm.sum, clump.enrich[perm==0], by=c("clinal.th", "gwas.th", "seas.th", "test"))

#calculate enrichment scores
perm.sum[,TT.enrich:=(TT-med.TT)/med.TT]
perm.sum[,prop.enrich:=(prop.test-med.prop.test)/med.prop.test]
perm.sum[,conc.enrich:=(conc.actual-med.prop.conc)/med.prop.conc]

perm.sum[prop.test>q.95.prop, enrich.prop:="greater"]
perm.sum[conc.actual>q.95.conc, enrich.conc:="greater"]
perm.sum[,gwas.th:=as.factor(gwas.th)]

#make plots


clump.enrich2<-fread("/scratch/pae3g/evolution/bergland2018_clineV2_enrich_intersectfirst_clump.txt")

#since intersections might be diffrent across permutations, need to calculate the fraction of sites that are clinal/seasonal rather than just the number


clump.enrich2[,prop.test:=TT/(TT+TF)]

perm.sum2<-clump.enrich2[perm!=0, .(med.TT=median(TT, na.rm=T),
                                    med.prop.test=median(prop.test, na.rm=T),
                                    med.prop.conc=median(conc.actual, na.rm=T),
                                    q.95.TT=quantile(TT, .95, na.rm=T),
                                    q.95.prop=quantile(prop.test, .95, na.rm=T),
                                    q.95.conc=quantile(conc.actual, .95, na.rm=T)),
                         .(clinal.th, gwas.th, seas.th, test)]

perm.sum2<-merge(perm.sum2, clump.enrich2[perm==0], by=c("clinal.th", "gwas.th", "seas.th", "test"))

#calculate enrichment scores
perm.sum2[,TT.enrich:=(TT-med.TT)/med.TT]
perm.sum2[,prop.enrich:=(prop.test-med.prop.test)/med.prop.test]
perm.sum2[,conc.enrich:=(conc.actual-med.prop.conc)/med.prop.conc]

perm.sum2[prop.test>q.95.prop, enrich.prop:="greater"]
perm.sum2[conc.actual>q.95.conc, enrich.conc:="greater"]
perm.sum2[,gwas.th:=as.factor(gwas.th)]

perm.sum2[,gwas.th:=as.factor(gwas.th)]


a<-ggplot(perm.sum[test=="cline"&!is.na(prop.enrich)], aes(x=as.factor(gwas.th), y=-1*clinal.th, fill=prop.enrich))+
    geom_tile()+
    scale_fill_viridis(option="viridis", limits=c(-1,4))+
    geom_tile(data=perm.sum[test=="cline"&prop.enrich==Inf], aes(x=as.factor(gwas.th), y=-1*clinal.th), fill="grey85")+
    geom_point(data=perm.sum[test=="cline"&enrich.prop=="greater"], aes(x=as.factor(gwas.th), y=-1*clinal.th), color="orange")+
    labs(x="", y=expression("-log"[10]*"(clinal q)"), fill="enrichment")+
    scale_x_discrete(labels=c("FDR < 0.025", "FDR < 0.05", "FDR < 0.1", "FDR < 0.15", "FDR < 0.2"))+
    theme(axis.text.x=element_blank())+
    theme(legend.position="bottom")+
    lims(y=c(0.7,5))


b<-ggplot(perm.sum[test=="seasonal"&!is.na(prop.enrich)], aes(x=as.factor(gwas.th), y=-1*seas.th, fill=prop.enrich))+
    geom_tile()+
    scale_fill_viridis(option="viridis", limits=c(-1,4))+
    geom_tile(data=perm.sum[test=="seasonal"&prop.enrich==Inf], aes(x=as.factor(gwas.th), y=-1*seas.th), fill="grey85")+
    geom_point(data=perm.sum[test=="seasonal"&enrich.prop=="greater"], aes(x=as.factor(gwas.th), y=-1*seas.th), color="orange")+
    labs(x="", y=expression("-log"[10]*"(seasonal q)"), fill="enrichment")+
    scale_x_discrete(labels=c("FDR < 0.025", "FDR < 0.05", "FDR < 0.1", "FDR < 0.15", "FDR < 0.2"))+
    theme(legend.position="bottom")+
    theme(axis.text.x=element_blank())+
    lims(y=c(0.7,5))

d<-get_legend(a)
title <- ggdraw() + draw_label("Bergland 2014", fontface='bold')



e<-ggplot(perm.sum2[test=="cline"&!is.na(prop.enrich)], aes(x=as.factor(gwas.th), y=-1*clinal.th, fill=prop.enrich))+
    geom_tile()+
    scale_fill_viridis(option="viridis", limits=c(-1,4))+
    geom_tile(data=perm.sum2[test=="cline"&prop.enrich==Inf], aes(x=as.factor(gwas.th), y=-1*clinal.th), fill="grey85")+
    geom_point(data=perm.sum2[test=="cline"&enrich.prop=="greater"], aes(x=as.factor(gwas.th), y=-1*clinal.th), color="orange")+
    labs(x="", y="", fill="enrichment")+
    scale_x_discrete(labels=c("FDR < 0.025", "FDR < 0.05", "FDR < 0.1", "FDR < 0.15", "FDR < 0.2"))+
    theme(axis.text.x=element_blank())+
    theme(axis.title.y=element_blank())+
    theme(legend.position="bottom")+
    lims(y=c(0.7,5))


f<-ggplot(perm.sum2[test=="seasonal"&!is.na(prop.enrich)], aes(x=as.factor(gwas.th), y=-1*seas.th, fill=prop.enrich))+
    geom_tile()+
    scale_fill_viridis(option="viridis", limits=c(-1,4))+
    geom_tile(data=perm.sum2[test=="seasonal"&prop.enrich==Inf], aes(x=as.factor(gwas.th), y=-1*seas.th), fill="grey85")+
    geom_point(data=perm.sum2[test=="seasonal"&enrich.prop=="greater"], aes(x=as.factor(gwas.th), y=-1*seas.th), color="orange")+
    labs(x="", y="", fill="enrichment")+
    scale_x_discrete(labels=c("FDR < 0.025", "FDR < 0.05", "FDR < 0.1", "FDR < 0.15", "FDR < 0.2"))+
    theme(legend.position="bottom")+
    theme(axis.text.x=element_blank())+
    theme(axis.title.y=element_blank())+
    lims(y=c(0.7,5))


h<-get_legend(e)
title <- ggdraw() + draw_label("Bergland 2014", fontface='bold')

title2 <- ggdraw() + draw_label("Machado 2019", fontface='bold')

#now plot polygenic scores

c.14<-fread("/scratch/pae3g/evolution/index_snps_PRS_cline2014.txt")
c.14[,avg.sum:=sum.Z.Z/n]
c.14[,rank:=frank(avg.sum), .(top)]
c.14[perm==0][order(rank)]

c.19<-fread("/scratch/pae3g/evolution/index_snps_PRS_cline2019.txt")
c.19[,avg.sum:=sum.Z.Z/n]
c.19[,rank:=frank(avg.sum), .(top)]
c.19[perm==0][order(rank)]


s.19<-fread("/scratch/pae3g/evolution/index_snps_PRS_seas2019.txt")
s.19[,avg.sum:=sum.Z.Z/n]
s.19[,rank:=frank(avg.sum), .(top)]
s.19[perm==0][order(rank)]

s.14<-fread("/scratch/pae3g/evolution/index_snps_PRS_seas2014.txt")
s.14[,avg.sum:=sum.Z.Z/n]
s.14[,rank:=frank(avg.sum), .(top)]
s.14[perm==0][order(rank)]


a2<-ggplot(c.14[perm!=0], aes(x=as.factor(top), y=avg.sum))+
    geom_quasirandom(method="smiley", color="grey80", size=0.5)+
    geom_point(data=c.14[perm==0], aes(x=as.factor(top), y=avg.sum), color="lightseagreen", size=1)+
    scale_x_discrete(labels=c("FDR < 0.025", "FDR < 0.05", "FDR < 0.1", "FDR < 0.15", "FDR < 0.2"))+
    labs(x="", y="clinal polygenic score")+theme(axis.text.x=element_blank())


b2<-ggplot(c.19[perm!=0], aes(x=as.factor(top), y=avg.sum))+
    geom_quasirandom(method="smiley", color="grey80", size=0.5)+
    geom_point(data=c.19[perm==0], aes(x=as.factor(top), y=avg.sum), color="lightseagreen", size=1)+
    scale_x_discrete(labels=c("FDR < 0.025", "FDR < 0.05", "FDR < 0.1", "FDR < 0.15", "FDR < 0.2"))+
    labs(x="", y="")+theme(axis.text.x=element_blank())




c2<-ggplot(s.14[perm!=0], aes(x=as.factor(top), y=avg.sum))+
    geom_quasirandom(method="smiley", color="grey80", size=0.5)+
    geom_point(data=s.14[perm==0], aes(x=as.factor(top), y=avg.sum), color="lightseagreen", size=1)+
    scale_x_discrete(labels=c("FDR < 0.025", "FDR < 0.05", "FDR < 0.1", "FDR < 0.15", "FDR < 0.2"))+
    labs(x="Index SNPs", y="seasonal polygenic score")+
    theme(axis.text.x=element_text(angle=45, hjust=1))

d2<-ggplot(s.19[perm!=0], aes(x=as.factor(top), y=avg.sum))+
    geom_quasirandom(method="smiley", color="grey80", size=0.5)+
    geom_point(data=s.19[perm==0], aes(x=as.factor(top), y=avg.sum), color="lightseagreen", size=1)+
    scale_x_discrete(labels=c("FDR < 0.025", "FDR < 0.05", "FDR < 0.1", "FDR < 0.15", "FDR < 0.2"))+
    labs(x="Index SNPs", y="")+
    theme(axis.text.x=element_text(angle=45, hjust=1))+theme(axis.title.y=element_blank())


right<-plot_grid(title2,
                 e+theme(legend.position="none"),
                 f+theme(legend.position="none"),
                 h+theme(legend.position="none"), 
                 b2,
                 d2,
                 ncol=1, rel_heights=c(0.03,0.2, 0.2,0.05, 0.2, 0.3), align="v", axis="l")
left<-plot_grid(title,
                a+theme(legend.position="none"),
                b+theme(legend.position="none"),
                d,
                a2,
                c2,
                ncol=1, rel_heights=c(0.03,0.2, 0.2,0.05, 0.2, 0.3), align="v", axis="l", labels=c("", "A", "B","" ,"C", "D"))


pdf("/scratch/pae3g/Figure5_final.pdf", height=8.75, width=6)

plot_grid(left, right, nrow=1,  rel_widths=c(0.49, 0.51), align="h", axis="b")
dev.off()



###############################
#### FIGURE 6 #################
###############################

library(data.table)
library(cowplot)
library(SNPRelate)
library(foreach)

seas<-fread("/scratch/pae3g/evolution/seas_glm_switch_clean.txt")
gwas<-fread(paste0("/scratch/pae3g/evolution/gwas_p_score_inv_id_perm", 0, ".txt"))

gwas<-merge(gwas, seas, by=c("chr", "pos"))

gwas[,rank.s:=frank(seas.p)/length(gwas$seas.p)]
gwas[,rank.g:=frank(gwas.p)/length(gwas$gwas.p)]

#plot of quantile values for seasonal and gwas

qgvsqs<-ggplot(gwas, aes(x=-log10(rank.g), y=-log10(rank.s)) )+geom_point(color="grey40", size=0.5)+
    geom_point(data=gwas[pos==9234730],aes(x=-log10(rank.g), y=-log10(rank.s)),size=2, color="lightseagreen")+labs(x=expression("-log"[10]*"(GWAS q)"), y=expression("-log"[10]*"(seas. q)"))+ theme(legend.position = "none")

#save jpeg of qgvsqs to add to pdf lataer

jpeg("/scratch/pae3g/Figure7_qgvqs.jpeg", height=2.66, width=3.85, res=1200, units="in")
qgvsqs
dev.off()


#plot effect of this snp

geno <- snpgdsOpen(filename = "/scratch/pae3g/genome-reconstruction/final3.vcf.gds", allow.fork=T)
phenos <- fread("/scratch/pae3g/phenos/phenos_062018.txt")
#pull genotypes
snpgenos.all <- snpgdsGetGeno(geno, snp.id = 926448,with.id=TRUE)

snpgenos<-data.table(snpgenos.all$genotype)
snpgenos[,sample.id:=snpgenos.all$sample.id]
dat=merge(phenos, snpgenos, by="sample.id")

a<-snpgdsSNPList(gdsobj=geno)

info<-data.table(snp.id=a$snp.id,
                 chr=a$chromosome,
                 pos=a$pos)

binomial_smooth <- function(...) {
    geom_smooth(method = "glm", method.args = list(family = "binomial"), se=FALSE)
}

effect<-ggplot(data=dat[!is.na(V1)], aes(x=temp.rack.cal, y=diapause.bin9, color=as.factor(V1)))+
    geom_jitter(data=dat[!is.na(V1)], aes(x=temp.rack.cal, y=diapause.bin9, color=as.factor(V1)), height = 0.05, alpha=0.5, size=0.2)+labs(x="Temperature °C", y="Probability\nof diapause", color="")+
    binomial_smooth()+scale_y_continuous(breaks=c(0,1))+
    theme(legend.position = c(0.5, 0.7))+
    scale_color_manual(values=c("#F8766D", "#7CAE00", "#00BFC4"), labels=c("fall/fall", "fall/spring", "spring/spring"))

snpgdsClose(geno)

#FST in recombinant outbred populations (workstation)

#open merged vcf

geno<-snpgdsOpen("/scratch/pae3g/evolution/inbred_and_dgrp.vcf.gds")

samps<-snpgdsSummary(geno)[["sample.id"]]

a<-snpgdsSNPList(gdsobj=geno)


info<-data.table(snp.id=a$snp.id,
                 chr=a$chromosome,
                 pos=a$pos)

#name the lines in two  recombinant outbred  populations
lines=c("12BME10-154", "12LN10-70", "Mayaguana_Mayaguana_43_18", "Ithica_I23", "12BME10-111", "12BME10-247", "12LN10-35", "12LN10-38", "12BME10-246", "12LN10-84", "Ithica_I26", "line_105", "line_109", "line_301", "line_313", "line_382", "line_129", "line_177", "line_181", "line_189")
d<-c("12BME10-270", "12LN10-95", "Freeport_GrandBahamasWest_33_11", "Ithica_I17", "12BME10-218", "12BME10-229", "12LN10-25", "12LN10-13", "12BME10-130", "12LN10-22", "Ithica_I22", "line_101", "line_142", "line_555", "line_646", "line_818", "line_375", "line_437", "line_712", "line_373")

group=as.factor(rep(c("A", "B"), each=length(lines)))

#calculate Fst between groups
c<-snpgdsFst(geno, population=group, sample.id=c(lines, d), with.id=T, autosome.only=F)
c.tab<-data.table(snp.id=c$snp.id,
                  fst=c$FstSNP)

#merge and assign colors for plot
c.tab<-merge(c.tab, info, by="snp.id")
c.tab[,color:=1]
c.tab[chr=="2R"|chr=="3R",color:=2]


#manhattan this so need a scale to translate snpid to chromosome
c.tab<-c.tab[chr!=4]
c.tab[,plot.id:=1:nrow(c.tab)]
chr.table<-c.tab[chr!="4",.(max.id=max(plot.id)), .(chr)]

chr.table<-foreach(i=1:5, .combine="rbind")%do%{
    temp<-chr.table[i]
    temp[,midpoint:=(chr.table[i,max.id]+chr.table[i-1, max.id])/2]
    return(temp)
}
chr.table[chr=="2L", midpoint:=max.id/2]
chr.table[,fdr.group:=NA]

#mir snp is not in this dataset so need to find what its snp id would be
info[chr=="2R"&pos>9234730]

#fst plot

fst<-ggplot(c.tab, aes(x=plot.id, y=fst, color=as.factor(color)))+ 
    geom_point(size=0.5)+
    labs(x=NULL, y=expression("F"["ST"]))+
    scale_color_manual(values=c("black", "grey"))+
    scale_x_continuous(breaks=chr.table$midpoint, labels=c("2L", "2R", "3L", "3R", "X"))+
    geom_point(aes(x=676410, y=1), color="lightseagreen", size=2) +
    theme(legend.position="none")


#try saving fst as a small jpeg to add into a pdf of other parts
jpeg("/scratch/pae3g/Figure7_fst.jpg", height=2.64, width=2.7, res=1200, units="in")
fst
dev.off()

### seasonal signal of mir snp####
library(data.table)
library(cowplot)

#load in processed data from this SNP from Machado et al 2018

mir.melt<-fread("/scratch/pae3g/evolution/mirsnpfreqs.csv")
season<-ggplot(mir.melt[set=="Core20"], aes( x=season, y=1-value, group=pop_name, color=season))+
    geom_point()+
    scale_color_manual(values=c("#F8766D","#00BFC4"), guide=F)+
    geom_line(data=mir.melt[set=="Core20"&switch==F], aes( x=season, y=1-value, group=pop_name), color="grey")+
    geom_line(data=mir.melt[set=="Core20"&switch==T], aes( x=season, y=1-value, group=pop_name), color="grey", linetype="dashed")+
    labs(x="", y="Allele Frequency") 


### quantitative complementaiton experiment
qdat<-fread("/scratch/pae3g/evolution/quant_comp_merged.txt")


#summary stats

dat.sum<-qdat[,.(prop.d9=sum(diapause.bin9)/.N, n=.N), .(cohort, genotype, Geno, season)]

dat.sum[,se.prop.d9:=sqrt((prop.d9*(1-prop.d9))/n)]

qc<-ggplot(dat.sum, aes(x=Geno, y=prop.d9, color=season, group=season))+geom_point()+geom_errorbar(data=dat.sum, aes(x=Geno, ymin=prop.d9-se.prop.d9, ymax=prop.d9+se.prop.d9), width=0.25, size=.5)+geom_line(size=.5)+labs(x="genetic background", y="Prop. diapause", color=NULL)+lims(y=c(.6, 1))+theme(axis.text.x=element_text(angle=45, hjust=1))+ theme(legend.position = c(0.05, 0.9))


#recombinant outbred populations

y<-fread("/scratch/pae3g/evolution/mir184_rop_052719.csv")
y[ovary>=10|egg>0, diapause.bin9:=0]
y[ovary<10&egg==0, diapause.bin9:=1]

y[ovary>=8|egg>0, diapause.bin:=0]
y[ovary<8&egg==0, diapause.bin:=1]
y[,vial:=as.factor(vial)]


y.sum<-y[,.(prop.d9=sum(diapause.bin9, na.rm=T)/.N, 
            prop.d=sum(diapause.bin, na.rm=T)/.N,
            n=.N), .(pop)]


y.sum[,se.prop.d9:=sqrt((prop.d9*(1-prop.d9))/n)]

rop<-ggplot(y.sum, aes(x=pop, y=prop.d9, color=pop))+geom_point()+geom_errorbar(data=y.sum, aes(x=pop, ymin=prop.d9-se.prop.d9, ymax=prop.d9+se.prop.d9, color=pop), width=0.5, size=.5)+lims(y=c(.6,1))+labs(y="Prop. diapause", x="")+scale_x_discrete(labels=c("fall", "spring"))+guides(color=F)

library(lubridate)
m<-fread("/scratch/pae3g/evolution/mir184_functional_022118.csv")
m[eggs==0&ovary<10, diapause.bin9:=1]
m[eggs>0|ovary>=10, diapause.bin9:=0]
m[short.geno=="c", geno:="Δ/CyO"]
m[short.geno=="w", geno:="Δ/Δ"]
m[,start:=mdy(eclosion)]
m[,stop:=mdy(freeze)]
m[,age:=stop-start]

#only look at batch 2 as it was better controlled
m.sum<-m[batch==2,.(prop.d9=sum(diapause.bin9)/.N, 
                    med.eggs=as.double(median(eggs)),
                    c=1/quantile(eggs, .75, na.rm=T),
                    mad.eggs=mad(eggs, constant = 1/quantile(eggs, .75, na.rm=T)),
                    n=.N),
         .(geno, temp, age)]
m.sum[,se.prop.d9:=sqrt(prop.d9*(1-prop.d9)/n)]



del.diap<-ggplot(m.sum, aes(x=as.factor(temp), y=prop.d9, color=geno, group=geno))+geom_point(position=position_dodge(width=0.5))+geom_errorbar(data=m.sum, aes(x=as.factor(temp), ymin=prop.d9-se.prop.d9, ymax=prop.d9+se.prop.d9, color=geno),position=position_dodge(width=0.5), width=0.5)+labs(x="temperature (°C)", y="Prop. diapause", color=NULL)+scale_colour_manual(values=c("#7CAE00", "#C77CFF"))+theme(axis.text.x=element_text(angle=45, hjust=1))+theme(legend.position = c(0.35, 0.9))


del.eggs<-ggplot(m.sum, aes(x=as.factor(temp), y=med.eggs, fill=geno))+geom_bar(stat="identity", position="dodge")+labs(x="temperature (°C)", y="Med. # eggs", fill=NULL )+scale_fill_manual(values=c("#7CAE00", "#C77CFF"))+theme(axis.text.x=element_text(angle=45, hjust=1))



blank<-ggplot()
##another version with mutant stuff
row1<-plot_grid(blank, effect, nrow=1, labels=c("A", "B"), align="h")
row2<-plot_grid(season,blank, rop, nrow=1, labels=c("C", "D", "E"), align="h", rel_widths=c(.35, .35, 0.3))
row3<-plot_grid(qc, 
                del.diap, 
                del.eggs+theme(legend.position="none"), 
                nrow=1, labels=c("F", "G", "H"), rel_widths = c(0.35, 0.35, 0.3), align="h")

pdf("/scratch/pae3g/Figure7_alt_parat.pdf", height=8, width=7.75)
plot_grid(row1, row2, row3, nrow=3, align="v", axis="l")
dev.off()


summary(glm(diapause.bin9~geno, data=m[batch==2&temp==25], family="binomial"))
summary(glm(diapause.bin9~geno, data=m[batch==2&temp==18], family="binomial"))

summary(glm(diapause.bin9~geno, data=m[batch==2&temp==12], family="binomial"))
summary(glm(diapause.bin9~geno, data=m[batch==2&temp==10.8], family="binomial"))


wilcox.test(eggs~geno, data=m[batch==2&temp==25])
wilcox.test(eggs~geno, data=m[batch==2&temp==18])#****
wilcox.test(eggs~geno, data=m[batch==2&temp==12])#*****
wilcox.test(eggs~geno, data=m[batch==2&temp==10.8])




#####################################
####### FIGURE 7 ####################
#####################################
library(data.table)
library(cowplot)
library(lme4)

f0<-fread("~/Box Sync/manuscripts/dryad_upload/mesocosm_F0.csv")
f0[,collection.date:=mdy(Date)]
f0[Eggs==0&Ovary<10, diapause.bin9:=1]
f0[Eggs>0|(Eggs==0&Ovary>=10), diapause.bin9:=0]
f0<-f0[!is.na(collection.date)&!is.na(Cage)]
f0[Cage%in%c(1:6), food:="fruit"]
f0[Cage%in%c(7:11), food:="cornmeal-mollasses"]

f0.sum<-f0[,.(prop.d9=sum(diapause.bin9, na.rm=T)/.N, 
              med.eggs=as.double(median(Eggs, na.rm=T)), 
              q75Eggs=quantile(Eggs, .75, na.rm=T),
              mad.eggs=mad(x=Eggs, constant=1/quantile(Eggs,.75, na.rm=T)),
              n=.N), .(collection.date, food)]
#add SE
f0.sum[,se.prop.d9:=sqrt(prop.d9*(1-prop.d9)/n)]


a<-ggplot(f0.sum[collection.date!="2019-01-25"&food=="fruit"], aes(x=collection.date, y=prop.d9, color="fruit-reared\noutdoor F0s")) +
    geom_point()+
    geom_line()+
    labs(x=NULL, y="Proportion of flies with\ndiapause-like ovaries", color="", linetype="") +
    geom_errorbar(data=f0.sum[collection.date!="2019-01-25"&food=="fruit"], aes(x=collection.date, ymin=prop.d9-se.prop.d9, ymax=prop.d9+se.prop.d9, color="fruit-reared\noutdoor F0s"), width=0.1)+
    geom_point(data=f0.sum[food=="cornmeal-mollasses"], aes(x=collection.date, y=prop.d9, color = "molasses-reared\noutdoor F0s"))+
    geom_errorbar(data=f0.sum[food=="cornmeal-mollasses"], aes(x=collection.date, ymin=prop.d9-se.prop.d9, ymax=prop.d9+se.prop.d9, color = "molasses-reared\noutdoor F0s"), width=0.1)+
    #geom_hline(aes(yintercept=dat.sum.date[collection.date=="2019-01-25", prop.d9], color="lab-reared indoor F0s", linetype="dashed"))+
    geom_hline(aes(yintercept = f0.sum[collection.date=="2019-01-25", prop.d9], linetype="indoor-reared F0s"), linetype="dashed", color="grey")+
    scale_color_manual(values=c("forestgreen", "saddlebrown"))+
    scale_x_date(date_breaks = "1 month", date_labels="%b", limits = as.Date(c('2018-06-15','2018-12-15')))

#stats for figure legend
summary(glmer(diapause.bin9~as.factor(collection.date)+(1|Cage), data=f0[collection.date!="2019-01-25"&collection.date!="2018-10-15"], family='binomial'))

summary(glmer(diapause.bin9~as.factor(collection.date)+(1|Cage), data=f0[collection.date=="2018-10-16"|collection.date=="2018-10-15"], family='binomial'))



f2<-fread("~/Box Sync/manuscripts/dryad_upload/mesocosm_F2.csv")
f2.sum=f2[,.(prop.d9=sum(diapause.bin9, na.rm=T)/.N,
             med.eggs=as.double(median(as.integer(eggs), na.rm=T)),
             mad.eggs=mad(x=eggs, constant=1/quantile(eggs,.75, na.rm=T)),
             n=.N), 
          .(box,location,date, temp)]
f2.sum[,se.prop.d9:=sqrt(prop.d9*(1-prop.d9)/n)]

b<-ggplot()+ geom_point(data=f2.sum, aes(x=date, y=prop.d9, color=location, group=interaction(location,temp)))+
    geom_line(data=f2.sum, aes(x=date, y=prop.d9, color=location, linetype=temp, group=interaction(location, temp)))+
    geom_errorbar(data=f2.sum, aes(x=date, ymin=prop.d9-se.prop.d9, ymax=prop.d9+se.prop.d9, color = location), width=0.1)+ labs(x=NULL, y="Proportion of flies\nin diapause", color="Location", linetype="Temperature") +lims(y=c(0,1))+
    scale_color_manual(values=c("grey","forestgreen" ))+
    scale_x_date(date_breaks = "1 month", date_labels="%b", limits = as.Date(c('2018-06-15','2018-12-15')))

#stats for figure legend
library(lme4)
summary(glmer(diapause.bin9~as.factor(date)+(1|cage), data=f2[location=="indoor"&box==46], family='binomial'))
summary(glmer(diapause.bin9~as.factor(date)+(1|cage), data=f2[location=="outdoor"&box==46], family='binomial'))
summary(glmer(diapause.bin9~as.factor(date)+(1|cage), data=f2[location=="indoor"&box==15], family='binomial'))
summary(glmer(diapause.bin9~as.factor(date)+(1|cage), data=f2[location=="outdoor"&box==15], family='binomial'))


w<-fread("~/Dropbox/Helen-Priscilla/weather_data.csv")
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


pdf("~/Desktop/Figure_7.pdf", height=8, width=8)

plot_grid(a, b, tempplot, align = "v", nrow = 3, rel_heights = c(.3, .5, .2), labels=c("A", "B", "C"))

dev.off()




#############################
###### FIGURE 8 #############
#############################


library(data.table)
library(cowplot)
library(ggbeeswarm)

ihs<-fread("/mnt/pricey_2/priscilla/ihs_dgrp_clump200kb.txt")

ihs.sum<-ihs[,
             .(n=.N,
               min.ihs=min(iHS, na.rm=T), 
               med.ihs=median(iHS, na.rm=T), 
               max.ihs=max(iHS, na.rm=T)), 
             .(top, perm)]

ihs.melt<-melt(ihs.sum, id.vars=c("top", "perm", "n"))


dpgp<-fread("/mnt/pricey_2/priscilla/dpgp/pro_diapause_alleles_africa_clump200kb.txt")

#calculate the proprotion of SNPs present above given allele frequencies
dpgp.sum<-dpgp[,.(n=.N,
                  prop.01=sum(pro.diapause>0.01, na.rm=T)/.N,
                  prop.05=sum(pro.diapause>0.05, na.rm=T)/.N,
                  prop.1=sum(pro.diapause>0.1, na.rm=T)/.N,
                  prop.2=sum(pro.diapause>0.2, na.rm=T)/.N,
                  med.p=median(pro.diapause, na.rm=T),
                  max.p=max(pro.diapause, na.rm=T),
                  min.p=min(pro.diapause, na.rm=T)),
               .(top, perm)]
#melt this data table
dpgp.melt<-melt(dpgp.sum, id.vars=c("top", "perm"), measure.vars=c("prop.01", "prop.05", "prop.1", "prop.2", "med.p"))

#calculate .05 quantile of permutations for one-tailed test (prediction is deficit of SNPs)
dpgp.q<-dpgp.sum[perm!=0, .(q.05.prop.01=quantile(prop.01, .05),
                            q.05.prop.05=quantile(prop.05, .05),
                            q.05.prop.1=quantile(prop.1, .05),
                            q.05.prop.2=quantile(prop.2, .05),
                            q.05.med.p=quantile(med.p, .05)), 
                 .(top)]    
#merge quantiles with actual data
dpgp.q<-merge(dpgp.sum[perm==0], dpgp.q, by="top")

#see which ones pass
dpgp.q[prop.01<=q.05.prop.01]
dpgp.q[prop.05<=q.05.prop.05]
dpgp.q[prop.1<=q.05.prop.1]
dpgp.q[prop.2<=q.05.prop.2]
dpgp.q[med.p<=q.05.med.p]



#tajima's d
tj<-fread("/mnt/pricey_2/priscilla/tajimasd_clumps200kb.txt")

tj.sum<-tj[,.(med=median(as.numeric(tajimas.d), na.rm=T),
              max=max(as.numeric(tajimas.d), na.rm=T),
              min=min(as.numeric(tajimas.d), na.rm=T)), 
           .(perm, top)]

tj.q<-tj.sum[perm!=0,.(q95.max=quantile(max, .95, na.rm=T)),.(top)]
tj.sum.melt<-melt(tj.sum, id.vars=c("perm", "top"))

tj.sum.melt[,variable:=factor(tj.sum.melt$variable, levels=c("min", "med", "max"))]


ihs.plot<-ggplot(ihs.melt[perm!=0], aes(x=as.factor(top), y=value, color=variable))+
    geom_quasirandom(dodge.width = .8, method="smiley", size=0.5)+
    geom_point(data=ihs.melt[perm==0], aes(x=as.factor(top), y=value, group=variable), color="black" ,position=position_dodge(width=0.8))+
    labs(x="", y="iHS in DGRP", color="")+
    scale_color_discrete(labels=c("minimum", "median", "maximum"))+
    theme(axis.text.x=element_blank())


tj.plot<-ggplot(tj.sum.melt[perm!=0], aes(x=as.factor(top), y=value, color=variable))+
    geom_quasirandom(dodge.width = .8, method="smiley", size=0.5)+
    geom_point(data=tj.sum.melt[perm==0], aes(x=as.factor(top), y=value, group=variable), color="black" ,position=position_dodge(width=0.8))+
    labs(x="Index SNPs", y="Tajima's D", color="")+
    scale_color_discrete(labels=c("minimum", "median", "maximum"))+
    scale_x_discrete(labels=c("FDR < 0.025", "FDR < 0.05", "FDR < 0.1", "FDR < 0.15", "FDR < 0.2"))+
    theme(axis.text.x=element_text(angle=45, hjust=1))



dpgp.plot<-ggplot(dpgp.melt[perm!=0], aes(x=as.factor(top), y=value, color=variable))+
    geom_quasirandom(dodge.width = .8, method="smiley", size=0.5,)+
    geom_point(data=dpgp.melt[perm==0], aes(x=as.factor(top), y=value, group=variable), color="black" ,position=position_dodge(width=0.8))+
    labs(x="", y="Proportion or\nallele frequency", color="Zambia")+
    scale_color_discrete(labels=c("p > 0.01", "p > 0.05", "p > 0.1", "p > 0.2", "median p"))+
    theme(axis.text.x=element_blank())


pdf("/mnt/pricey_2/priscilla/Figure8.pdf", height=7, width=7)
plot_grid( dpgp.plot,ihs.plot, tj.plot, nrow=3, align="v", labels=c("A", "B", "C"), rel_heights = c(0.3, 0.3, 0.4))
dev.off()



###################################
#### FIGURE S1 ####################
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
#### FIGURE S2 ##################
#################################


#supplemental figure showing generation and swarm with each phenotype continues from the data table generated in Figure 2

#melt phenotypes into long format
p.melt<-melt(p, id.vars=c("id", "temp.rack.cal", "photoperiod", "generation", "swarm"), measure.vars=c("diapause.bin9", "diapause.bin", "eggP"), variable.name="pheno", value.name="diapause")
p.melt[,pheno:=factor(p.melt$pheno, levels=c("eggP", "diapause.bin9", "diapause.bin"))]

p.melt[pheno=="eggP", stage:="Stage 14"]
p.melt[pheno=="diapause.bin", stage:="Stage 8"]
p.melt[pheno=="diapause.bin9", stage:="Stage 10"]

p.melt[,stage:=factor(p.melt$stage, levels=c("Stage 8", "Stage 10", "Stage 14"))]

#Generation plot
a.1<-ggplot(p.melt, aes(x=temp.rack.cal, y=diapause, color=as.factor(generation)))+
    binomial_smooth()+
    geom_point(data=p.melt, aes(x=temp.rack.cal, y=diapause),color="grey", alpha=0.1)+
    labs(color="Generation", x="Temperature °C", y="Probability of diapause")+facet_grid(.~stage)

#Cage plot
b.1<-ggplot(p.melt, aes(x=temp.rack.cal, y=diapause, color=as.factor(swarm)))+
    binomial_smooth()+
    geom_point(data=p.melt, aes(x=temp.rack.cal, y=diapause), color="grey", alpha=0.1)+
    labs(color="Population", x="Temperature °C", y="Probability of diapause")+facet_grid(.~stage)
#Final plot: effect of generation and cage on phenotypes

pdf("~/Box Sync/manuscripts/FigureS2.pdf", height=6, width=8)
plot_grid(a.1, b.1, nrow=2, labels=c("A", "B"))
dev.off()




#################################
#### FIGURE S3 ######## #########
#################################

library(data.table)
library(cowplot)


data<-fread("~/Box Sync/hybridSwarm/yeast_supplementation.csv")

p.melt<-melt(data, id.vars=c( "temp.rack.cal", "box", "treatment", "cage"), measure.vars=c("diapause9", "diapause", "eggp"), variable.name="pheno", value.name="diapause")

p.melt[pheno=="eggp", stage:="Stage 14"]
p.melt[pheno=="diapause", stage:="Stage 8"]
p.melt[pheno=="diapause9", stage:="Stage 10"]

p.melt[,stage:=factor(p.melt$stage, levels=c("Stage 8", "Stage 10", "Stage 14"))]

pdf("~/Box Sync/manuscripts/FigureS3.pdf", heigh=3.5, width=8)
ggplot(p.melt, aes(x=temp.rack.cal,y=diapause, color=treatment, linetype=cage))+
        geom_jitter(height=0.1,alpha=0.5)+
        binomial_smooth()+
        labs(x="Temperature (°C)",y="Probability of diapause")+
        scale_y_continuous(breaks=c(0,1))+
    labs(color=NULL)+
    facet_wrap(~stage)
dev.off()

summary(aov(glm(diapause~temp.rack.cal+cage+yeast, family=binomial, data=data)))
summary(aov(glm(diapause9~temp.rack.cal+cage+yeast, family=binomial, data=data)))
summary(aov(glm(eggp~temp.rack.cal+cage+yeast, family=binomial, data=data)))




#########################
## Figure S4 ############
#########################

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

#################################
### FIGURE S5 ###################
#################################


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


######################
### FIGURE S6 ########  
######################


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
##### FIGURE S7 ##############
##############################

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


######################################
### FIGURE S8 ########################
######################################

library(data.table)
library(cowplot)
pcaOut=fread("/scratch/pae3g/evolution/final3.vcf.PCA.csv")

kary.ag<-fread("/scratch/pae3g/evolution/final3.vcf.karytypecalls.csv", header=T)
kary.ag.melt=melt(kary.ag, measure.vars=c("chr2L", "chr2R", "chr3L", "chr3R"), id.vars="sample.id")
names(kary.ag.melt)=c("sample.id", "chrs", "kary")

pcaOut=merge(pcaOut, kary.ag.melt, by=c("sample.id", "chrs"), all=TRUE)
pcaOut[is.na(kary), kary:="std"]
pcaOut[chrs=="chr2L", chr:="2L"]
pcaOut[chrs=="chr3L", chr:="3L"]
pcaOut[chrs=="chr2R", chr:="2R"]
pcaOut[chrs=="chr3R", chr:="3R"]
pcaOut[chrs=="chrX", chr:="X"]
pcaOut[chrs=="all", chr:="all"]
pcaOut[,kary2:=ifelse(kary=="std", "Standard", "Inversion")]
pcaOut[,chr:=factor(chr, levels=c("all", "2L", "2R", "3L", "3R", "X"))]


#reorder karyotype levels
pcaOut[,kary:=factor(kary, levels=c("std", "std;In(2L)t" , "std;In(2R)Ns" , "std;In(3R)C", "std;In(3R)Mo" , "std;In(3R)Payne","In(2L)t","In(2R)Ns","In(3R)C", "In(3R)Mo","In(3R)Payne", "In(3R)Payne;In(3R)C", "In(3R)Payne;In(3R)Mo" ))]



pdf("/scratch/pae3g/pc_kary.pdf", height=6, width=10)
ggplot(pcaOut, aes(x=pc1, y=pc2, color=actual.cage, shape=kary2))+geom_point()+facet_wrap(~chr, scales="free")+labs(color="cage", shape="karyotype")
dev.off()



########################
##### FIGURE S9 ########
########################

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

pdf("/scratch/pae3g/IBS_GRM.pdf", height=8, width=8)
ggplot(a.melt, aes(x=id, y=Var2, fill=value))+
    geom_tile()+
    scale_fill_viridis(option="viridis")+ 
    #scale_fill_gradientn(colours = terrain.colors(20))+
    theme(line = element_blank(),
          axis.text = element_blank(),
          title = element_blank())+
    labs(fill="IBS")

dev.off()

################################################
############# FIGURE S10 #######################
################################################

library(data.table)
library(cowplot)

#averaged data
g<-fread("/mnt/sammas_storage/bergland-lab/Priscilla/perm_summaries/adaptive_perms_perm0_pvalues_simple.txt")
#single imputation
a<-fread("/mnt/sammas_storage/bergland-lab/Priscilla/gwas_imputation_permutations/genesis_diapause.bin9~generation+temp.rack.cal+photoperiod_both_draw29_perm0_replaced.txt")

#plots
gd<-ggplot(a, aes(x=Score.pval))+geom_histogram(bins=1000)+labs(x="Single imputation GENESIS P-values")
gp<-ggplot(g, aes(x=avg.pval))+geom_histogram(bins=1000)+labs(x="GENESIS P-values averaged across imputations")
ep<-ggplot(g, aes(x=pval.percentile))+geom_histogram(bins=1000)+labs(x="Empirical P-values")
ge<-ggplot(g, aes(x=-1*log10(avg.pval), y=-1*log10(pval.percentile), color=as.factor(maf.bin)))+geom_point(size=0.5)+geom_abline(intercept=0, slope=1, linetype="dashed")+labs(x=expression("-log"[10]*"(Avg. GENESIS P)"), y=expression("-log"[10]*"(Empirical P)"), color="MAF bin")

jpeg("/mnt/pricey_2/priscilla/Pvalue_dist.jpg", height=8, width=6,res=1200, units="in")
plot_grid(gd, gp, ep, ge, ncol=1, labels=c("A", "B", "C", "D"), rel_heights = c(0.2, 0.2, .2, 0.4), align="v", axis="l")

dev.off()


###################################
########## FIGURE S11 #############
###################################

#continues from data table produced in Figure 3
p2L<-ggplot(y[perm!=0], aes(x=l.2L))+geom_histogram(binwidth=0.05, fill="grey85")+geom_vline(xintercept=y[perm==0,l.2L], color="lightseagreen", size=2)+labs(x="λ GC", title="Chr. 2L")+scale_x_continuous(limits=c(.5, 2.5), breaks=c(0.5,1, 1.5, 2,2.5))
p2R<-ggplot(y[perm!=0], aes(x=l.2R))+geom_histogram(binwidth=0.05, fill="grey85")+geom_vline(xintercept=y[perm==0,l.2R], color="lightseagreen", size=2)+labs(x="λ GC", title="Chr. 2R")+scale_x_continuous(limits=c(.5, 2.5), breaks=c(0.5,1, 1.5, 2,2.5))
p3L<-ggplot(y[perm!=0], aes(x=l.3L))+geom_histogram(binwidth=0.05, fill="grey85")+geom_vline(xintercept=y[perm==0,l.3L], color="lightseagreen", size=2)+labs(x="λ GC", title="Chr. 3L")+scale_x_continuous(limits=c(.5, 2.5), breaks=c(0.5,1, 1.5, 2,2.5))
p3R<-ggplot(y[perm!=0], aes(x=l.3R))+geom_histogram(binwidth=0.05, fill="grey85")+geom_vline(xintercept=y[perm==0,l.3R], color="lightseagreen", size=2)+labs(x="λ GC", title="Chr. 3R")+scale_x_continuous(limits=c(.5, 2.5), breaks=c(0.5,1, 1.5, 2,2.5))
pX<-ggplot(y[perm!=0], aes(x=l.X))+geom_histogram(binwidth=0.05, fill="grey85")+geom_vline(xintercept=y[perm==0,l.X], color="lightseagreen", size=2)+labs(x="λ GC", title="Chr. X")+scale_x_continuous(limits=c(.5, 2.5), breaks=c(0.5,1, 1.5, 2,2.5))



pdf("/scratch/pae3g/Figure_S12.pdf", height=2.5, width=12)
plot_grid(p2L, p2R, p3L, p3R, pX, nrow=1)
dev.off()




#########################
## FIGURE S12 ##########
#########################

library(data.table)
library(SNPRelate)
library(cowplot)

tim.geno<-fread("/scratch/pae3g/tim_genotypes.csv")

phenos<-fread("/nv/vol186/bergland-lab/Priscilla/phenos_062018.txt")
phenos[photoperiod==9, pp:="9L:15D"]
phenos[photoperiod==11, pp:="11L:13D"]
phenos[photoperiod==13, pp:="13L:11D"]
phenos[photoperiod==15, pp:="15L:9D"]
phenos[,pp:=factor(phenos$pp, levels=c("9L:15D", "11L:13D", "13L:11D", "15L:9D"))]
tim.genos<-merge(phenos, tim.geno, by="id")


#open up hybrid file to pull cpo genotypes
geno<-snpgdsOpen("/scratch/pae3g/genome-reconstruction/final3.vcf.gds")

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

left.top<-top<-ggplot(data=tim.genos[!is.na(final.geno)], aes(x=temp.rack.cal, y=diapause.bin9, color=as.factor(final.geno)))+geom_jitter(data=tim.genos[!is.na(final.geno)], aes(x=temp.rack.cal, y=diapause.bin9, color=as.factor(final.geno)), height = 0.1)+labs(x="Temperature °C", y="Diapause", color="tim genotype")+binomial_smooth()+scale_y_continuous(breaks=c(0,1))+scale_color_discrete(labels=c("s-tim/s-tim", "s-tim/ls-tim", "ls-tim/ls-tim"))

left.bottom<-ggplot(data=cpo[!is.na(cpo.geno)], aes(x=temp.rack.cal, y=diapause.bin9, color=as.factor(cpo.geno)))+geom_jitter(data=cpo[!is.na(cpo.geno)], aes(x=temp.rack.cal, y=diapause.bin9, color=as.factor(cpo.geno)), height = 0.1)+labs(x="Temperature °C", y="Diapause", color="cpo genotype")+binomial_smooth()+scale_y_continuous(breaks=c(0,1))+scale_color_discrete(labels=c("A/A", "A/T", "T/T"))


top<-ggplot(data=tim.genos[!is.na(final.geno)], aes(x=temp.rack.cal, y=diapause.bin9, color=as.factor(final.geno)))+geom_jitter(data=tim.genos[!is.na(final.geno)], aes(x=temp.rack.cal, y=diapause.bin9, color=as.factor(final.geno)), height = 0.1)+labs(x="Temperature °C", y="Diapause", color="tim genotype")+binomial_smooth()+scale_y_continuous(breaks=c(0,1))+facet_grid(.~pp)+scale_color_discrete(labels=c("s-tim/s-tim", "s-tim/ls-tim", "ls-tim/ls-tim"))


bottom<-ggplot(data=cpo[!is.na(cpo.geno)], aes(x=temp.rack.cal, y=diapause.bin9, color=as.factor(cpo.geno)))+geom_jitter(data=cpo[!is.na(cpo.geno)], aes(x=temp.rack.cal, y=diapause.bin9, color=as.factor(cpo.geno)), height = 0.1)+labs(x="Temperature °C", y="Diapause", color="cpo genotype")+binomial_smooth()+scale_y_continuous(breaks=c(0,1))+facet_grid(.~pp)+scale_color_discrete(labels=c("A/A", "A/T", "T/T"))

top.legend<-get_legend(left.top)
bottom.legend<-get_legend(left.bottom)
top.row<-plot_grid(left.top+theme(legend.position="none"),
                   top+theme(legend.position="none"),
                   top.legend,
                   labels=c("A", "B", ""),
                   rel_widths=c(.25, .6, .15),
                   nrow=1 )


bottom.row<-plot_grid(left.bottom+theme(legend.position="none"),
                   bottom+theme(legend.position="none"),
                   bottom.legend,
                   labels=c("C", "D", ""),
                   rel_widths=c(.25, .6, .15),
                   nrow=1 )

pdf("/scratch/pae3g/cpo_tim.pdf", height=6, width=10)
plot_grid(top.row, bottom.row, nrow=2)
dev.off()


#stats for figure legend
summary(glm(diapause.bin9~final.geno+temp.rack.cal+photoperiod+swarm.y+generation, data=tim.genos, family=binomial))
summary(glm(diapause.bin9~final.geno+temp.rack.cal+swarm.y+generation, data=tim.genos[photoperiod==9], family=binomial))
summary(glm(diapause.bin9~final.geno+temp.rack.cal+swarm.y+generation, data=tim.genos[photoperiod==11], family=binomial))
summary(glm(diapause.bin9~final.geno+temp.rack.cal+swarm.y+generation, data=tim.genos[photoperiod==13], family=binomial))
summary(glm(diapause.bin9~final.geno+temp.rack.cal+swarm.y+generation, data=tim.genos[photoperiod==15], family=binomial))

summary(glm(diapause.bin9~cpo.geno+temp.rack.cal+photoperiod+swarm+generation, data=cpo, family=binomial))
summary(glm(diapause.bin9~cpo.geno+temp.rack.cal+swarm+generation, data=cpo[photoperiod==9], family=binomial))
summary(glm(diapause.bin9~cpo.geno+temp.rack.cal+swarm+generation, data=cpo[photoperiod==11], family=binomial))
summary(glm(diapause.bin9~cpo.geno+temp.rack.cal+swarm+generation, data=cpo[photoperiod==13], family=binomial))
summary(glm(diapause.bin9~cpo.geno+temp.rack.cal+swarm+generation, data=cpo[photoperiod==15], family=binomial))



#############################
### FIGURE S13 ##############
#############################


library(data.table)
library(foreach)
library(cowplot)
library(viridis)
library(SNPRelate)
geno <- snpgdsOpen("/scratch/pae3g/genome-reconstruction/final2_draw1_replaced.vcf.gds", allow.fork=T)
a<-snpgdsSNPList(gdsobj=geno)

info<-data.table(snp.id=a$snp.id,
                 chr=a$chromosome,
                 pos=a$pos)

perm.snps<-foreach(perm=c(0,101:119), .errorhandling="remove")%do%{
    print(perm)
    gwas<-fread(paste0("/scratch/pae3g/evolution/gwas_p_score_inv_id_perm", perm, ".txt"))
    gwas<-merge(gwas, info, by=c("chr", "pos"))
    snps<-(gwas[order(gwas.p)][1:101,snp.id])
    ld<-snpgdsLDMat(geno, slide=-1,verbose=F,snp.id =snps , method='composite')$LD
    m<-data.table(melt(ld))
    m[,perm:=perm]
    return(m)
}

perm.snps<-rbindlist(perm.snps)

pdf("/scratch/pae3g/Figure_S13.pdf", height=7, width=8)
ggplot(perm.snps, aes(x=Var1, y=Var2, fill=value*value))+
    geom_tile()+
    scale_fill_viridis(option="viridis")+ 
    theme(line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())+
    labs( fill=expression("R"^2))+
    facet_wrap(~perm, scales="free")
dev.off()

################################################
########### FIGURE S14 #########################
################################################



library(data.table)
library(cowplot)
library(ggbeeswarm)
library(foreach)

clumps<-fread("/mnt/pricey_2/priscilla/index_snp_perms.txt")
pdf("/mnt/pricey_2/priscilla/index_snp_counts.pdf",height=6, width=6)

ggplot(clumps[perm!=0], aes(x=as.factor(top), y=log10(n)))+
    geom_quasirandom(method="smiley", color="grey40", size=0.5)+geom_point(data=clumps[perm==0], aes(x=as.factor(top), y=log10(n)), color="lightseagreen", size=2)+
    scale_x_discrete(labels=c("FDR < .025", "FDR < 0.05", "FDR < 0.1", "FDR < 0.15", "FDR < 0.2"))+
    labs(y=expression("log"[10]*"(# index SNPs)"), x="")
dev.off()




##########################
#### FIGURE S15 ##########
##########################

#continues from code to make Figure 5!

perm.sum[,excess:=TT-med.TT]
perm.sum2[,excess:=TT-med.TT]

#deal with 0s in excess
perm.sum[,log10excess:=sign(excess)*log10(abs(excess))]
perm.sum2[,log10excess:=sign(excess)*log10(abs(excess))]

perm.sum[is.na(log10excess), log10excess:=0]
perm.sum2[is.na(log10excess), log10excess:=0]

a<-ggplot(perm.sum[test=="cline"&TT>0], 
          aes(x=gwas.th, 
              y=-1*clinal.th, 
              fill=log10excess))+
    geom_tile()+
    scale_fill_viridis(option="viridis",limits=c(-.5, 3))+
    labs(x="", 
         y=expression("-log"[10]*"(clinal q)"), 
         fill=expression("log"[10]*"(#)"))+
    geom_point(data=perm.sum[test=="cline"&enrich.prop=="greater"], 
               aes(x=gwas.th, y=-1*clinal.th), 
               color="orange", 
               fill="orange")+
    theme(legend.position="bottom")+
    theme(axis.text.x=element_blank())+
    lims( y=c(0.7,5.1))

b<-ggplot(perm.sum2[test=="cline"&TT>0], 
          aes(x=gwas.th, 
              y=-1*clinal.th, 
              fill=log10excess))+
    geom_tile()+
    scale_fill_viridis(option="viridis",limits=c(-.5, 3))+
    labs(x="", 
         y="", 
         fill=expression("log"[10]*"(#)"))+
    geom_point(data=perm.sum2[test=="cline"&enrich.prop=="greater"], 
               aes(x=gwas.th, y=-1*clinal.th), 
               color="orange", 
               fill="orange")+
    theme(legend.position="bottom")+
    theme(axis.text.x=element_blank())+
    lims( y=c(0.7,5.1))


c<-ggplot(perm.sum[test=="seasonal"&TT>0], 
          aes(x=gwas.th, 
              y=-1*seas.th, 
              fill=log10excess))+
    geom_tile()+
    scale_fill_viridis(option="viridis",limits=c(-.5, 3))+
    labs(x="", 
         y=expression("-log"[10]*"(seasonal q)"), 
         fill=expression("log"[10]*"(#)"))+
    geom_point(data=perm.sum[test=="seasonal"&enrich.prop=="greater"], 
               aes(x=gwas.th, y=-1*seas.th), 
               color="orange", 
               fill="orange")+
    theme(legend.position="bottom")+
    lims( y=c(0.7,5.1))+
    scale_x_discrete(labels=c("FDR < 0.025", "FDR < 0.05", "FDR < 0.1", "FDR < 0.15", "FDR < 0.2"))+
    theme(axis.text.x=element_text(angle=45, hjust=1))


d<-ggplot(perm.sum2[test=="seasonal"&TT>0], 
          aes(x=gwas.th, 
              y=-1*seas.th, 
              fill=log10excess))+
    geom_tile()+
    scale_fill_viridis(option="viridis",limits=c(-.5, 3))+
    labs(x="", 
         y="", 
         fill=expression("log"[10]*"(#)"))+
    geom_point(data=perm.sum2[test=="seasonal"&enrich.prop=="greater"], 
               aes(x=gwas.th, y=-1*seas.th), 
               color="orange", 
               fill="orange")+
    theme(legend.position="bottom")+
    lims( y=c(0.7,5.1))+
    scale_x_discrete(labels=c("FDR < 0.025", "FDR < 0.05", "FDR < 0.1", "FDR < 0.15", "FDR < 0.2"))+
    theme(axis.text.x=element_text(angle=45, hjust=1))



h<-get_legend(a)
title <- ggdraw() + draw_label("Bergland 2014", fontface='bold')

title2 <- ggdraw() + draw_label("Machado 2019", fontface='bold')

left<-plot_grid(title, a+theme(legend.position="none"), c+theme(legend.position="none"), h, ncol=1, rel_heights=c(0.1, 0.4, 0.5, 0.1), labels=c("", "A", "B"), align="v")
right<-plot_grid(title2,b +theme(legend.position="none"), d+theme(legend.position="none"), ggplot(), ncol=1, rel_heights=c(0.1, 0.4, 0.5, 0.1), align="v")

pdf("/scratch/pae3g/FigureS15_excess.pdf", height=6, width=6)
plot_grid(left, right, ncol=2, align="h")
dev.off()



###########################################
### FIGURE S16#############################
###########################################

library(data.table)
library(cowplot)
library(viridis)

clump.enrich<-fread("/scratch/pae3g/evolution/index_snps_enrich_by_pop_seas2019.txt")

clump.enrich[,prop.test:=TT/(TT+TF)]

#summarize all permutations except 0, which is original ordering of data
perm.sum<-clump.enrich[perm!=0, .(med.TT=median(TT, na.rm=T),
                                  med.prop.test=median(prop.test, na.rm=T),
                                  q.95.TT=quantile(TT, .95, na.rm=T),
                                  q.95.prop=quantile(prop.test, .95, na.rm=T)),
                       .( gwas.th, seas.th, test, pop)]

#merge with actual data
perm.sum<-merge(perm.sum, clump.enrich[perm==0], by=c( "gwas.th", "seas.th", "pop", "test"))

#calculate enrichment scores
perm.sum[,TT.enrich:=(TT-med.TT)/med.TT]
perm.sum[,prop.enrich:=(prop.test-med.prop.test)/med.prop.test]

perm.sum[prop.test>q.95.prop, enrich.prop:="greater"]
perm.sum[,gwas.th:=as.factor(gwas.th)]

ids=fread("/scratch/pae3g/evolution/core20pops.csv")
perm.sum<-merge(perm.sum, ids, by="pop")

perm.sum<-merge()
#for clarity, set enrichment above 5 to 5
perm.sum[,prop.enrich.adj:=prop.enrich]
perm.sum[prop.enrich>5&prop.enrich<Inf, prop.enrich.adj:=5]

pdf("/scratch/pae3g/enrich_bypop_supp.pdf", height=8, width=7)
ggplot(perm.sum[!is.na(prop.enrich)], aes(x=as.factor(gwas.th), y=-1*seas.th, fill=prop.enrich.adj))+
    geom_tile()+
    scale_fill_viridis(option="viridis", limits=c(-1,5))+
    geom_tile(data=perm.sum[prop.enrich==Inf], aes(x=as.factor(gwas.th), y=-1*seas.th), fill="grey85")+
    geom_point(data=perm.sum[enrich.prop=="greater"], aes(x=as.factor(gwas.th), y=-1*seas.th), color="orange")+
    labs(x="", y=expression("-log"[10]*"(seas q)"), fill="enrichment")+
    scale_x_discrete(labels=c("FDR < 0.025", "FDR < 0.05", "FDR < 0.1", "FDR < 0.15", "FDR < 0.2"))+
    theme(axis.text.x=element_text(angle=45, hjust=1))+
    theme(legend.position="bottom")+
    #lims(y=c(0.7,5))
    facet_wrap(~full_name)
dev.off()



################################
##### FIGURE S17 ###############
################################

library(data.table)
library(foreach)
library(cowplot)
library(ggbeeswarm)

c<-fread("/scratch/pae3g/evolution/index_snps_PRS_bypop_seas2019.txt")

#switch signs of scores for switch populations
c[pop%in%c("TKA_14", "BHM_14", "LMA_14", "rd_12"), sum.logit:=-1*sum.logit]
c[,avg.logit:=sum.logit/n]
c[,rank.avg.logit:=frank(avg.logit), .(top, pop)]

c[perm==0&rank.avg.logit>91]
c[perm==0&rank.avg.logit>96]

pdf("/scratch/pae3g/FigureS17.pdf", height=4, width=6)
ggplot(c[perm!=0&pop%in%c("LMA_14", "TKA_14", "AGA_14")], aes(x=as.factor(top), y=avg.logit))+
    geom_quasirandom(method="smiley", color="grey80", size=0.5)+
    geom_point(data=c[perm==0&pop%in%c("LMA_14", "TKA_14", "AGA_14")], aes(x=as.factor(top), y=avg.logit), color="lightseagreen", size=1)+
    scale_x_discrete(labels=c("FDR < 0.025", "FDR < 0.05", "FDR < 0.1", "FDR < 0.15", "FDR < 0.2"))+
    theme(axis.text.x=element_text(angle=45, hjust=1))+
    labs(x="", y="seasonal polygenic score")+
    facet_wrap(~pop)
dev.off()
