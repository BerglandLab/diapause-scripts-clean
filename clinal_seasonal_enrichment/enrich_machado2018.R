#!/usr/bin/env Rscript


library(data.table)
library(foreach)
library(doMC)
registerDoMC(10)



enrich<-function(gwas, snpset, cline, seas, af, clinal.th.list, seas.th.list, seas.fixed.th, gwas.th.list, cline_polar, seas_polar, min.maf=0, max.maf=.5, output){
    print("generating datasets")
    #load data. gwas will just be the clumped snps
    gwas.dat<-fread(gwas)
    snpset<-fread(snps)
    #format plink output
    setnames(snpset, c("CHR", "BP"), c("chr", "pos"))
    snpset[,chr:=as.character(chr)]
    snpset[chr=="23", chr:="X"]
    #merge with gwas
 
    #read cline and season and calculate q-values
    cline.dat<-fread(cline)
    cline.dat<-merge(gwas.dat[,.(chr,pos)], cline.dat, by=c("chr", "pos"))
    seas.dat<-fread(seas)
    seas.dat<-merge(gwas.dat[,.(chr,pos)], seas.dat, by=c("chr", "pos"))
    b<-merge(cline.dat, seas.dat, by=c("chr", "pos"))
    b[,q.c:=log10(frank(clinal.p)/(.N+1))]
    b[,q.s:=log10(frank(seas.p)/(.N+1))]
    af.dat<-fread(af)
    gwas.dat<-merge(gwas.dat, snpset[,.(chr, pos)], by=c("chr", "pos"))
    #merge things
   
    b.g<-merge(gwas.dat, b, by=c("chr", "pos"))
    #merge in allele frequency data and calculate
    b.g<-merge(b.g, af.dat, by=c("chr", "pos"))
    b.g[,maf:=pmin(f.hat, 1-f.hat)]
    #filter based on minor allele frequency
    b.g<-b.g[maf>=min.maf&maf<=max.maf]
    #give ids to final dataset
    b.g[,id:=c(1:dim(b.g)[1])]
    setkey(b.g, chr, pos)
    #correct coefficients to make concordant. in alan's data, positive means that ALT allele is higher in north and fall. in my data, postiive mean REF allele increases diapause
    #change coefficient for latitude to make predicated signs all go in the same direction. BE SURE THAT THIS IS RIGHT!!!
    b.g[,clinal.beta:=cline_polar*clinal.beta]
    b.g[,seas.beta:=seas_polar*seas.beta]
    
    #loop through different gwas thresholds 
    
    #do clinal calculations if required. produce two tables that will get rbinded at the end
    if(!is.null(clinal.th.list)){
        #loop through clinal thersholds
        clinal<-foreach(clinal.q.thresh = clinal.th.list, .errorhandling="remove")%dopar%{
            #permform clinal enrichment/concordance tets
            cline.enrich<-data.table(TT=sum(b.g$q.c<=clinal.q.thresh, na.rm=T),
                                     TF=sum(b.g$q.c>clinal.q.thresh, na.rm=T),
                                     TT.conc=sum(b.g$q.c<=clinal.q.thresh&sign(b.g$gwas.Score)==sign(b.g$clinal.beta), na.rm=T),
                                     TF.conc=sum(b.g$q.c>clinal.q.thresh&sign(b.g$gwas.Score)==sign(b.g$clinal.beta), na.rm=T),
                                     clinal.th=clinal.q.thresh,
                                     gwas.th=top,
                                     seas.th=1,
                                     test="cline")
            
            #make clinal and seasonal table if needed at fixed seasonal threshold
            if(!is.null(clinal.th.list)&!is.null(seas.th.list)){
                cline.seas.enrich<-data.table(TT=sum(b.g$q.c<=clinal.q.thresh  & b.g$q.s<seas.fixed.th, na.rm=T),
                                              TF=sum(b.g$q.c>clinal.q.thresh  & b.g$q.s<seas.fixed.th, na.rm=T),
                                              TT.conc=sum(b.g$q.c<=clinal.q.thresh & b.g$q.s<seas.fixed.th&sign(b.g$gwas.Score)==sign(b.g$clinal.beta)&sign(b.g$gwas.Score)==sign(b.g$seas.beta), na.rm=T),
                                              TF.conc=sum(b.g$q.c>clinal.q.thresh & b.g$q.s<seas.fixed.th&sign(b.g$gwas.Score)==sign(b.g$clinal.beta)&sign(b.g$gwas.Score)==sign(b.g$seas.beta), na.rm=T),
                                              clinal.th=clinal.q.thresh,
                                              gwas.th=top,
                                              seas.th=seas.fixed.th,
                                              test="cline+season")
                
            }
            return(rbind(cline.enrich, cline.seas.enrich))
        }
        
        clinal<-rbindlist(clinal)
    } else{
        clinal<-dat.table()
    }
    
    #separate analysis for seasonal 
    if(!is.null(seas.th.list)){ 
        print('testing seasonal enrichment')
        seasonal<-foreach(seas.q.thresh = seas.th.list, .errorhandling="remove")%dopar%{
            return(data.table(TT=sum(b.g$q.s<=seas.q.thresh, na.rm=T),
                              TF=sum(b.g$q.s>seas.q.thresh, na.rm=T),
                              TT.conc=sum(b.g$q.s<=seas.q.thresh&sign(b.g$gwas.Score)==sign(b.g$seas.beta), na.rm=T),
                              TF.conc=sum(b.g$q.s>seas.q.thresh&sign(b.g$gwas.Score)==sign(b.g$seas.beta), na.rm=T),
                              clinal.th=1,
                              gwas.th=top,
                              seas.th=seas.q.thresh,
                              test="seasonal"
            ))
        }
        seasonal<-rbindlist(seasonal)
    } else{
        seasonal<-data.table()
    }
    
gwas.tab<-rbind(clinal, seasonal)
gwas.tab[,conc.actual:=TT.conc/TT]
gwas.tab[,min.maf:=min.maf]
gwas.tab[,max.maf:=max.maf]
gwas.tab[,perm:=perm]
return(gwas.tab)
}

#use function to loop through perms and thresholds
clump.enrich<-foreach(perm =c(0,101:200)) %do%{
    print(perm)
    a<-foreach (top = c(101, 490, 5316, 25250, 67914)) %do% {
        
        cline<-"/scratch/pae3g/evolution/east_coast_cline_V2_clean.txt"
        seas<-"/scratch/pae3g/evolution/seas_glm_switch_clean.txt"
        freqs<-"/scratch/pae3g/evolution/east_coast_cline_V2_allele_freqs.txt"
        gwas<-paste0("/scratch/pae3g/evolution/gwas_p_score_inv_id_perm", perm, ".txt")
        snps<-paste0("/nv/vol186/bergland-lab/Priscilla/perm_summaries/perm", perm, "_top", top, "_ld0.4_200kb.clumped")
        
        enrich(gwas=gwas,
               snps=snps,
               cline=cline,
               seas=seas,
               af=freqs,
               clinal.th.list=seq(from=-5, to=-.8, by=0.1), 
               seas.th.list=seq(from=-5, to=-.8, by=0.1), 
               seas.fixed.th=log10(.05), 
               gwas.th.list=seq(from=-5, to=-.8, by=0.1),
               cline_polar=-1, 
               seas_polar=-1, 
               min.maf=0,
               max.maf=0.5,
               output=paste0("/scratch/pae3g/evolution/bergland2018_clineV2_all_alleles_perm", perm, "_top",top ))
    }
    return(rbindlist(a))
}

clump.enrich=rbindlist(clump.enrich)



write.table(clump.enrich, "/scratch/pae3g/evolution/bergland2018_clineV2_enrich_intersectfirst_clump.txt", quote=F, row.names=F, sep="\t")



###analyze and plot
library(data.table)
library(cowplot)
library(viridis)

clump.enrich<-fread("/scratch/pae3g/evolution/bergland2018_clineV2_enrich_intersectfirst_clump.txt")

#since intersections might be diffrent across permutations, need to calculate the fraction of sites that are clinal/seasonal rather than just the number


clump.enrich[,prop.test:=TT/(TT+TF)]

perm.sum<-clump.enrich[perm!=0, .(med.TT=median(TT, na.rm=T),
                                  med.prop.test=median(prop.test, na.rm=T),
                                  med.prop.conc=median(conc.actual, na.rm=T),
                                  q.025.TT=quantile(TT, .025, na.rm=T),
                                  q.975.TT=quantile(TT, .975, na.rm=T),
                                  q.025.prop=quantile(prop.test, .025, na.rm=T),
                                  q.975.prop=quantile(prop.test, .975, na.rm=T),
                                  q.025.conc=quantile(conc.actual, .025, na.rm=T),
                                  q.975.conc=quantile(conc.actual, .975, na.rm=T)),
                       .(clinal.th, gwas.th, seas.th, test)]

perm.sum<-merge(perm.sum, clump.enrich[perm==0], by=c("clinal.th", "gwas.th", "seas.th", "test"))

perm.sum[,TT.enrich:=(TT-med.TT)/med.TT]
perm.sum[,prop.enrich:=(prop.test-med.prop.test)/med.prop.test]
perm.sum[,conc.enrich:=(conc.actual-med.prop.conc)/med.prop.conc]

perm.sum[TT<=q.025.TT, enrich.TT:="less"]
perm.sum[TT>=q.975.TT, enrich.TT:="greater"]
perm.sum[prop.test<q.025.prop, enrich.prop:="less"]
perm.sum[prop.test>q.975.prop, enrich.prop:="greater"]
perm.sum[conc.actual<q.025.conc, enrich.conc:="less"]
perm.sum[conc.actual>q.975.conc, enrich.conc:="greater"]
perm.sum[,gwas.th:=as.factor(gwas.th)]



a<-ggplot(perm.sum[test=="cline"&!is.na(prop.enrich)], aes(x=as.factor(gwas.th), y=-1*clinal.th, fill=prop.enrich))+
    geom_tile()+
    scale_fill_viridis(option="viridis", limits=c(-1,4))+
    geom_tile(data=perm.sum[test=="cline"&prop.enrich==Inf], aes(x=as.factor(gwas.th), y=-1*clinal.th), fill="grey85")+
    geom_point(data=perm.sum[test=="cline"&enrich.prop=="greater"], aes(x=as.factor(gwas.th), y=-1*clinal.th), color="orange")+
    geom_point(data=perm.sum[test=="cline"&enrich.prop=="lesser"], aes(x=as.factor(gwas.th), y=-1*clinal.th), color="white")+
    labs(x="", y=expression("-log"[10]*"(Clinal q)"), fill="enrichment", title="clinal")+
    scale_x_discrete(labels=c("FDR < 0.025", "FDR < 0.05", "FDR < 0.1", "FDR < 0.15", "FDR < 0.2"))+
    theme(axis.text.x=element_text(angle=45, hjust=1))


b<-ggplot(perm.sum[test=="seasonal"&!is.na(prop.enrich)], aes(x=as.factor(gwas.th), y=-1*seas.th, fill=prop.enrich))+
    geom_tile()+
    scale_fill_viridis(option="viridis", limits=c(-1,4))+
    geom_tile(data=perm.sum[test=="seasonal"&prop.enrich==Inf], aes(x=as.factor(gwas.th), y=-1*seas.th), fill="grey85")+
    geom_point(data=perm.sum[test=="seasonal"&enrich.prop=="greater"], aes(x=as.factor(gwas.th), y=-1*seas.th), color="orange")+
    geom_point(data=perm.sum[test=="seasonal"&enrich.prop=="lesser"], aes(x=as.factor(gwas.th), y=-1*seas.th), color="white")+
    labs(x="", y=expression("-log"[10]*"(Seasonal q)"), fill="enrichment",title = "seasonal")+
    scale_x_discrete(labels=c("FDR < 0.025", "FDR < 0.05", "FDR < 0.1", "FDR < 0.15", "FDR < 0.2"))+
    theme(axis.text.x=element_text(angle=45, hjust=1))

c<-ggplot(perm.sum[test=="cline+season"&!is.na(prop.enrich)], aes(x=as.factor(gwas.th), y=-1*clinal.th, fill=prop.enrich))+
    geom_tile()+
    scale_fill_viridis(option="viridis", limits=c(-1,4))+
    geom_tile(data=perm.sum[test=="cline+season"&prop.enrich==Inf], aes(x=as.factor(gwas.th), y=-1*clinal.th), fill="grey85")+
    geom_point(data=perm.sum[test=="cline+season"&enrich.prop=="greater"], aes(x=as.factor(gwas.th), y=-1*clinal.th), color="orange")+
    geom_point(data=perm.sum[test=="cline+season"&enrich.prop=="lesser"], aes(x=as.factor(gwas.th), y=-1*clinal.th), color="white")+
    labs(x="", y=expression("-log"[10]*"(Clinal q)"), fill="enrichment", title="clinal; seasonal q < 0.05")+
    scale_x_discrete(drop=F, labels=c("FDR < 0.025", "FDR < 0.05", "FDR < 0.1", "FDR < 0.15", "FDR < 0.2"))+
    theme(axis.text.x=element_text(angle=45, hjust=1))

d<-get_legend(a)

pdf("/scratch/pae3g/clinal_seasonal_clump_enrich.pdf", height=4, width=8)
plot_grid(a+theme(legend.position="none"),
          b+theme(legend.position="none"),
          c+theme(legend.position="none"),
          d, nrow=1, rel_widths=c(0.3, 0.3, 0.3, 0.1))
dev.off()


#concordance
d<-ggplot(perm.sum[test=="cline"], aes(x=as.factor(gwas.th), y=clinal.th, fill=conc.enrich))+
    geom_tile()+
    scale_fill_viridis(option="viridis")+
    geom_point(data=perm.sum[test=="cline"&enrich.conc=="greater"], aes(x=as.factor(gwas.th), y=clinal.th), color="orange")+
    geom_point(data=perm.sum[test=="cline"&enrich.conc=="lesser"], aes(x=as.factor(gwas.th), y=clinal.th), color="white")


e<-ggplot(perm.sum[test=="seasonal"], aes(x=as.factor(gwas.th), y=seas.th, fill=conc.enrich))+
    geom_tile()+
    scale_fill_viridis(option="viridis")+
    geom_point(data=perm.sum[test=="seasonal"&enrich.conc=="greater"], aes(x=as.factor(gwas.th), y=seas.th), color="orange")+
    geom_point(data=perm.sum[test=="seasonal"&enrich.conc=="lesser"], aes(x=as.factor(gwas.th), y=seas.th), color="white")


f<-ggplot(perm.sum[test=="cline+season"], aes(x=as.factor(gwas.th), y=clinal.th, fill=conc.enrich))+
    geom_tile()+
    scale_fill_viridis(option="viridis")+
    geom_point(data=perm.sum[test=="cline+season"&enrich.conc=="greater"], aes(x=as.factor(gwas.th), y=clinal.th), color="orange")+
    geom_point(data=perm.sum[test=="cline+season"&enrich.conc=="lesser"], aes(x=as.factor(gwas.th), y=clinal.th), color="white")


plot_grid(d,e,f, nrow=1)