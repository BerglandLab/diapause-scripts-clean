

library(data.table)
library(foreach)
library(cowplot)
library(doMC)
registerDoMC(20)

seas.th.list=seq(from=-5, to=-.8, by=0.1)

ids=fread("/scratch/pae3g/evolution/all_popinfo_formap_set.csv")
pops=unique(ids[set=="Core20",pop_name])

pop.test<-foreach(pop=pops)%do%{
    print(pop)
    clump.enrich<-foreach(perm =c(0,101:200)) %dopar% { #cycle through original data and 100 permutations
        print(perm)
        a<-foreach (top = c(101, 490, 5316, 25250, 67914)) %do% { #these numbers correspond to the top # of snps in our different FDR cutoffs
            
            gwas<-fread(paste0("/scratch/pae3g/evolution/gwas_p_score_inv_id_perm", perm, ".txt"))
            snps<-fread(paste0("/nv/vol186/bergland-lab/Priscilla/perm_summaries/perm", perm, "_top", top, "_ld0.4_200kb.clumped"))
            seas<-fread(paste0("/scratch/pae3g/evolution/fisher_exactJ_", pop,".coef_minp2.txt"))
                        setnames(seas, "chrom", "chr")
                        setnames(snps, c("CHR", "BP"), c("chr", "pos"))
                        snps[,chr:=as.character(chr)]
                        snps[chr=="23", chr:="X"]
                        gwas<-merge(gwas, seas, by=c("chr", "pos"))
                        
                        
                        #quantile rank and normalize fisher P values in merged datset
                        gwas[,q.s:=log10(frank(minp2)/(length(minp2)+1))]
                        gwas<-merge(gwas, snps[, .(chr, pos)], by=c("chr", "pos"))
                        
                        #return sums
                        seasonal<-foreach(seas.q.thresh = seas.th.list, .errorhandling="remove")%dopar%{
                            return(data.table(TT=sum(gwas$q.s<=seas.q.thresh, na.rm=T),
                                              TF=sum(gwas$q.s>seas.q.thresh, na.rm=T),
                                              gwas.th=top,
                                              seas.th=seas.q.thresh,
                                              test="seasonal",
                                              pop=pop,
                                              perm=perm
                            ))
                        }
                        return(rbindlist(seasonal))
        }
        return(rbindlist(a))
    }
    
    clump.enrich<-rbindlist(clump.enrich)
    clump.enrich[,prop.test:=TT/(TT+TF)]
    
    return(clump.enrich)
}

pop.test<-rbindlist(pop.test)
write.table(pop.test, "/scratch/pae3g/evolution/index_snps_enrich_by_pop_seas2019.txt", quote=F, sep="\t", row.names=F)

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


#what is the really enriched snp in MA12?

perm=0
pop="MA_12"
top=101
gwas<-fread(paste0("/scratch/pae3g/evolution/gwas_p_score_inv_id_perm", perm, ".txt"))
snps<-fread(paste0("/nv/vol186/bergland-lab/Priscilla/perm_summaries/perm", perm, "_top", top, "_ld0.4_200kb.clumped"))
seas<-fread(paste0("/scratch/pae3g/evolution/fisher_exactJ_", pop,".coef_minp2.txt"))
setnames(seas, "chrom", "chr")
setnames(snps, c("CHR", "BP"), c("chr", "pos"))
snps[,chr:=as.character(chr)]
snps[chr=="23", chr:="X"]
gwas<-merge(gwas, seas, by=c("chr", "pos"))
gwas[,q.s:=log10(frank(minp2)/(length(minp2)+1))]
gwas<-merge(gwas, snps[, .(chr, pos)], by=c("chr", "pos"))

#chr3L:8205502, psoitive score so ref increases diapause

load(file="/scratch/pae3g/evolution/dat.Rdata")
setnames(dat, "pop", 'name')

d<-dat[pos==8205502]

ids=fread("/scratch/pae3g/evolution/all_popinfo_formap_set.csv")

d<-merge(d, ids, by="name")
d[,switch:=ifelse(pop_name%in%c("TKA_14", "BHM_14", "LMA_14", "rd_12"), T, F)]
d[switch==T&season=="fall", season.switch:="spring"]
d[switch==T&season=="spring", season.switch:="fall"]
d[switch==T, season:=season.switch]
d.melt<-melt(d[,.(pop_name, season, af,set)], id.vars=c("pop_name", "season", "set"), measure.vars=c("af"))
d.melt[,switch:=ifelse(pop_name%in%c("TKA_14", "BHM_14", "LMA_14", "rd_12"), T, F)]
d.melt[,season:=factor(d.melt$season, levels=c("spring", "fall"))]
# 


ggplot(d.melt[set=="Core20"], aes( x=season, y=1-value, group=pop_name, color=season))+
    geom_point()+
    scale_color_manual(values=c("#00BFC4","#F8766D"), guide=F)+
    geom_line(data=d.melt[set=="Core20"&switch==F], aes( x=season, y=1-value, group=pop_name), color="grey")+
    geom_line(data=d.melt[set=="Core20"&switch==T], aes( x=season, y=1-value, group=pop_name), color="grey", linetype="dashed")+
    labs(x="", y="Allele Frequency", title="3L:8205502 (Mass. '12)")

#goes strongly in wrong direction in 3 pops and mildly in right direction in most others


#what is the really enriched snp in MA14?

perm=0
pop="LMA_14"
top=101
gwas<-fread(paste0("/scratch/pae3g/evolution/gwas_p_score_inv_id_perm", perm, ".txt"))
snps<-fread(paste0("/nv/vol186/bergland-lab/Priscilla/perm_summaries/perm", perm, "_top", top, "_ld0.4_200kb.clumped"))
seas<-fread(paste0("/scratch/pae3g/evolution/fisher_exactJ_", pop,".coef_minp2.txt"))
setnames(seas, "chrom", "chr")
setnames(snps, c("CHR", "BP"), c("chr", "pos"))
snps[,chr:=as.character(chr)]
snps[chr=="23", chr:="X"]
gwas<-merge(gwas, seas, by=c("chr", "pos"))
gwas[,q.s:=log10(frank(minp2)/(length(minp2)+1))]
gwas<-merge(gwas, snps[, .(chr, pos)], by=c("chr", "pos"))
gwas[q.s<(-3)]
#chr2L:15883645, neg score so alt allele increase diapause

e<-dat[pos==15883645]


e<-merge(e, ids, by="name")
e[,switch:=ifelse(pop_name%in%c("TKA_14", "BHM_14", "LMA_14", "rd_12"), T, F)]
e[switch==T&season=="fall", season.switch:="spring"]
e[switch==T&season=="spring", season.switch:="fall"]
e[switch==T, season:=season.switch]
e.melt<-melt(e[,.(pop_name, season, af,set)], id.vars=c("pop_name", "season", "set"), measure.vars=c("af"))
e.melt[,switch:=ifelse(pop_name%in%c("TKA_14", "BHM_14", "LMA_14", "rd_12"), T, F)]
e.melt[,season:=factor(d.melt$season, levels=c("spring", "fall"))]
# 


ggplot(e.melt[set=="Core20"], aes( x=season, y=value, group=pop_name, color=season))+
    geom_point()+
    scale_color_manual(values=c("#00BFC4","#F8766D"), guide=F)+
    geom_line(data=e.melt[set=="Core20"&switch==F], aes( x=season, y=value, group=pop_name), color="grey")+
    geom_line(data=e.melt[set=="Core20"&switch==T], aes( x=season, y=value, group=pop_name), color="grey", linetype="dashed")+
    labs(x="", y="Pro-diapause allele frequency", title="2L:15883645 (Mass. '14)")





#what is the really enriched snp in PA_10?

perm=0
pop="PA_10"
top=490
gwas<-fread(paste0("/scratch/pae3g/evolution/gwas_p_score_inv_id_perm", perm, ".txt"))
snps<-fread(paste0("/nv/vol186/bergland-lab/Priscilla/perm_summaries/perm", perm, "_top", top, "_ld0.4_200kb.clumped"))
seas<-fread(paste0("/scratch/pae3g/evolution/fisher_exactJ_", pop,".coef_minp2.txt"))
setnames(seas, "chrom", "chr")
setnames(snps, c("CHR", "BP"), c("chr", "pos"))
snps[,chr:=as.character(chr)]
snps[chr=="23", chr:="X"]
gwas<-merge(gwas, seas, by=c("chr", "pos"))
gwas[,q.s:=log10(frank(minp2)/(length(minp2)+1))]
gwas<-merge(gwas, snps[, .(chr, pos)], by=c("chr", "pos"))
gwas[q.s<(-3)]
#chr2L:15883645, neg score so alt allele increase diapause

f<-dat[pos==12972118]


f<-merge(f, ids, by="name")
f[,switch:=ifelse(pop_name%in%c("TKA_14", "BHM_14", "LMA_14", "rd_12"), T, F)]
f[switch==T&season=="fall", season.switch:="spring"]
f[switch==T&season=="spring", season.switch:="fall"]
f[switch==T, season:=season.switch]
f.melt<-melt(f[,.(pop_name, season, af,set)], id.vars=c("pop_name", "season", "set"), measure.vars=c("af"))
f.melt[,switch:=ifelse(pop_name%in%c("TKA_14", "BHM_14", "LMA_14", "rd_12"), T, F)]
f.melt[,season:=factor(d.melt$season, levels=c("spring", "fall"))]
# 


ggplot(f.melt[set=="Core20"], aes( x=season, y=value, group=pop_name, color=season))+
    geom_point()+
    scale_color_manual(values=c("#00BFC4","#F8766D"), guide=F)+
    geom_line(data=f.melt[set=="Core20"&switch==F], aes( x=season, y=value, group=pop_name), color="grey")+
    geom_line(data=f.melt[set=="Core20"&switch==T], aes( x=season, y=value, group=pop_name), color="grey", linetype="dashed")+
    labs(x="", y="Pro-diapause allele frequency", title="2R:12972118 (Penn. '10)")

#goes strongly in wrong direction in 3 pops and mildly in right direction in most others
