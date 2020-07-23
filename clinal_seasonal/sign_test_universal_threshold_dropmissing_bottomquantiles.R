library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(foreach)
library(doMC)
registerDoMC(16)
library(viridis)
library(ggbeeswarm)
library(Rmisc)

print("bergland 2014")

load("/scratch/pae3g/oldscratch_recovered/evolution/6d_data.Rdata")
b<-as.data.table(p)
b[,maf:=pmin(f.hat, 1-f.hat)]
b<-b[maf>=0.05]
b<-b[clinal.beta<3&clinal.beta>(-3)]
b<-b[sfsfsfX.beta<3&sfsfsfX.beta>(-3)]
b[,clinal.q:=log10(frank(clinal.p)/(length(clinal.p)+1))]
b[,sfsfsfX.q:=log10(frank(sfsfsfX.p)/(length(sfsfsfX.p)+1))]

files<-fread("/scratch/pae3g/genome-reconstruction/universal_input.txt")



#read in files
y<-foreach(pop=files$V1[1:3000], phenotype=files$V2[1:3000],draw=files$V3[1:3000], perm=files$V4[1:3000], model=files$V6[1:3000])%dopar%{
    print(paste(draw, perm, sep=","))
    #read gwas adn do some fixes
    load(paste("/scratch/pae3g/revisions/genesis_", phenotype, "_draw", draw, "_perm", perm, "_pop", pop, "_" , model, "_allsnpgrm_wolbachia_dropmissing.Rdat", sep=""))
    gwas<-assoc.results
    gwas[,maf:=pmin(freq, 1-freq)]
    gwas<-gwas[maf>=0.05]
    gwas[,q:=log10(frank(Score.pval)/(length(Score.pval)+1))]
    gwas<-merge(gwas, b[,.(chr, pos, clinal.beta, clinal.p,clinal.q,  maf, sfsfsfX.p, sfsfsfX.beta, sfsfsfX.q)], by=c("chr", "pos"))
    gwas[,TT.cline:=ifelse(sign(Score)==1 & sign(clinal.beta)==(-1), T, F)]
    gwas[,TF.cline:=ifelse(sign(Score)==1 & sign(clinal.beta)==1, T, F)]
    gwas[,FT.cline:=ifelse(sign(Score)==(-1) & sign(clinal.beta)==(-1), T, F)]
    gwas[,FF.cline:=ifelse(sign(Score)==(-1) & sign(clinal.beta)==1, T, F)]

    gwas[,TT.seas:=ifelse(sign(Score)==1 & sign(sfsfsfX.beta)==1, T, F)]
    gwas[,TF.seas:=ifelse(sign(Score)==1 & sign(sfsfsfX.beta)==(-1), T, F)]
    gwas[,FT.seas:=ifelse(sign(Score)==(-1) & sign(sfsfsfX.beta)==1, T, F)]
    gwas[,FF.seas:=ifelse(sign(Score)==(-1) & sign(sfsfsfX.beta)==(-1), T, F)]

    #bergland cline is based on alternate alleles and needs to be flipped. 
    gwas[,ps.cline:=-1*Score.Stat*clinal.beta]
    #bergland positive clinal beta means alternate allele is higher in the fall (same direction as gwas)
    gwas[,ps.seas:=1*Score.Stat*sfsfsfX.beta]


    cline.top<-foreach (top = seq(from=-5, to=0, by=1)) %do% {
        a<-gwas[q>top,.(
            n=.N,
            TT=sum(TT.cline, na.rm=T),
            TF=sum(TF.cline, na.rm=T),
            FT=sum(FT.cline, na.rm=T),
            FF=sum(FF.cline, na.rm=T),
            poly=sum(ps.cline)
        )]
        a[,top:=top]
    }
    cline.top<-rbindlist(cline.top)
    cline.top[,or:=as.numeric(TT)*as.numeric(FF)/(as.numeric(FT)*as.numeric(TF))]
    cline.top[, perm:=perm]
    cline.top[,draw:=draw]
    cline.top[,test:="clinal"]
    cline.top[,pheno:=phenotype]
    cline.top[,model:=model]
    cline.top[,pop:=pop]

    

    seas.top<-foreach (top = seq(from=-5, to=0, by=1)) %do% {
        a<-gwas[q>top,.(
            n=.N,
            TT=sum(TT.seas, na.rm=T),
            TF=sum(TF.seas, na.rm=T),
            FT=sum(FT.seas, na.rm=T),
            FF=sum(FF.seas, na.rm=T),
            poly=sum(ps.seas)
        )]
        a[,top:=top]
    }

    seas.top<-rbindlist(seas.top)
    seas.top[,or:=as.numeric(TT)*as.numeric(FF)/(as.numeric(FT)*as.numeric(TF))]
    seas.top[, perm:=perm]
    seas.top[,draw:=draw]
    seas.top[,test:="seasonal"]
    seas.top[,pheno:=phenotype]
    seas.top[,model:=model]
    seas.top[,pop:=pop]
    

    return(rbind(cline.top, seas.top))
}

y<-rbindlist(y)

save(y, file="/scratch/pae3g/revisions/evolution/bergland2014_sign_universal_threshold_dropmissing_bottom_quantiles.Rdata")

# print("machado 2019")
# 
# cline<-fread("/scratch/pae3g/revisions/evolution/east_coast_cline_V2_clean.txt")
# seas<-fread("/scratch/pae3g/revisions/evolution/seas_glm_NOTswitch_clean.txt")   
# freqs<-fread("/scratch/pae3g/oldscratch_recovered/evolution/east_coast_cline_V2_allele_freqs.txt")
# 
# x<-merge(cline, seas, by=c("chr", "pos"))
# x<-merge(x, freqs, by=c("chr", "pos"))
# 
# x<-x[f.hat>=0.05&f.hat<=0.95]
# x<-x[clinal.beta<3&clinal.beta>(-3)]
# x<-x[seas.beta<3&seas.beta>(-3)]
# 
# 
# 
# 
# 
# y<-foreach(pop=files$V1[1:3000], phenotype=files$V2[1:3000],draw=files$V3[1:3000], perm=files$V4[1:3000], model=files$V6[1:3000], .errorhandling="remove")%dopar%{
#     print(paste(draw, perm, sep=","))
#     #read gwas adn do some fixes
#     load(paste("/scratch/pae3g/revisions/genesis_", phenotype, "_draw", draw, "_perm", perm, "_pop", pop, "_" , model, "_allsnpgrm_wolbachia_dropmissing.Rdat", sep=""))
#     gwas<-assoc.results
#     gwas[,maf:=pmin(freq, 1-freq)]
#     gwas<-gwas[maf>=0.05]
#     gwas[,q:=log10(frank(Score.pval)/(length(Score.pval)+1))]
#     gwas<-merge(gwas, x[,.(chr, pos, clinal.beta, clinal.p, seas.beta)], by=c("chr", "pos"))
#     gwas[,TT.cline:=ifelse(sign(Score)==1 & sign(clinal.beta)==(-1), T, F)]
#     gwas[,TF.cline:=ifelse(sign(Score)==1 & sign(clinal.beta)==1, T, F)]
#     gwas[,FT.cline:=ifelse(sign(Score)==(-1) & sign(clinal.beta)==(-1), T, F)]
#     gwas[,FF.cline:=ifelse(sign(Score)==(-1) & sign(clinal.beta)==1, T, F)]
#     
#     gwas[,TT.seas:=ifelse(sign(Score)==1 & sign(seas.beta)==(-1), T, F)]
#     gwas[,TF.seas:=ifelse(sign(Score)==1 & sign(seas.beta)==(1), T, F)]
#     gwas[,FT.seas:=ifelse(sign(Score)==(-1) & sign(seas.beta)==(-1), T, F)]
#     gwas[,FF.seas:=ifelse(sign(Score)==(-1) & sign(seas.beta)==(1), T, F)]
#     
#     gwas[,ps.cline:=-1*Score.Stat*clinal.beta]
#     #note that the seasonal polygenic score should have had a *-1 in it and needs to be flipped for plotting
#     gwas[,ps.seas:=1*Score.Stat*seas.beta]
#     
#     
#     cline.top<-foreach (top = seq(from=-5, to=0, by=1)) %do% {
#         a<-gwas[q<top,.(
#             n=.N,
#             TT=sum(TT.cline, na.rm=T),
#             TF=sum(TF.cline, na.rm=T),
#             FT=sum(FT.cline, na.rm=T),
#             FF=sum(FF.cline, na.rm=T),
#             poly=sum(ps.cline)
#         )]
#         a[,top:=top]
#     }
#     
#     cline.top<-rbindlist(cline.top)
#     cline.top[,or:=as.numeric(TT)*as.numeric(FF)/(as.numeric(FT)*as.numeric(TF))]
#     cline.top[, perm:=perm]
#     cline.top[,draw:=draw]
#     cline.top[,test:="clinal"]
#     cline.top[,pheno:=phenotype]
#     cline.top[,model:=model]
#     cline.top[,pop:=pop]
#     
#     seas.top<-foreach (top = seq(from=-5, to=0, by=1)) %do% {
#         a<-gwas[q<top,.(
#             n=.N,
#             TT=sum(TT.seas, na.rm=T),
#             TF=sum(TF.seas, na.rm=T),
#             FT=sum(FT.seas, na.rm=T),
#             FF=sum(FF.seas, na.rm=T),
#             poly=sum(ps.seas)
#         )]
#         a[,top:=top]
#     }
#     
#     seas.top<-rbindlist(seas.top)
#     seas.top[,or:=as.numeric(TT)*as.numeric(FF)/(as.numeric(FT)*as.numeric(TF))]
#     seas.top[, perm:=perm]
#     seas.top[,draw:=draw]
#     seas.top[,test:="seasonal"]
#     seas.top[,pheno:=phenotype]
#     seas.top[,model:=model]
#     seas.top[,pop:=pop]
#     
#     
#     return(rbind(cline.top, seas.top))
# }
# 
# y<-rbindlist(y)
# 
# save(y, file="/scratch/pae3g/revisions/evolution/bergland2019_sign_universal_threshold_dropmissing.Rdata")
# 
# 
# print("individual populations")
# load("/scratch/pae3g/oldscratch_recovered/evolution/core20delta.rdat")
# pops<-names(deltas)
# 
# 
# 
# 
# y<-foreach(pop=files$V1[1:3000], phenotype=files$V2[1:3000],draw=files$V3[1:3000], perm=files$V4[1:3000], model=files$V6[1:3000], .errorhandling="remove")%dopar%{
#     print(paste(draw, perm, sep=","))
#     #read gwas adn do some fixes
#     load(paste("/scratch/pae3g/revisions/genesis_", phenotype, "_draw", draw, "_perm", perm, "_pop", pop, "_" , model, "_allsnpgrm_wolbachia_dropmissing.Rdat", sep=""))
#     gwas<-assoc.results
#     gwas[,maf:=pmin(freq, 1-freq)]
#     gwas<-gwas[maf>=0.05]
#     gwas[,q:=log10(frank(Score.pval)/(length(Score.pval)+1))]
# 
#     pop.test<-foreach(p=pops)%do%{
#         print(p)
#         seas.top<-foreach (top = c(-5:0)) %do% {
#             
#             l<-merge(gwas[q<=top], deltas[[p]][is.finite(diff.logit)], by=c("chr", "pos"))
#             l[,TT:=ifelse(sign(Score.Stat)==1 & sign(diff.logit)==(1), T, F)]
#             l[,TF:=ifelse(sign(Score.Stat)==(-1) & sign(diff.logit)==1, T, F)]
#             l[,FT:=ifelse(sign(Score.Stat)==(1) & sign(diff.logit)==(-1), T, F)]
#             l[,FF:=ifelse(sign(Score.Stat)==(-1) & sign(diff.logit)==(-1), T, F)]
#             l[,ps:=Score.Stat*diff.logit]
#             l[,top:=top]
#             return(l[,.(or=sum(as.numeric(TT), na.rm=T)*sum(as.numeric(FF), na.rm=T)/(sum(as.numeric(TF), na.rm=T)*sum(as.numeric(FT), na.rm=T)), 
#                         poly=sum(ps)), .(population, top)])
#             
#         }
#         return(rbindlist(seas.top))
#     }
#     
#     
#     
#     pop.test<-rbindlist(pop.test)
#     pop.test[,pheno:=phenotype]
#     pop.test[,pop:=pop]
#     pop.test[, perm:=perm]
#     pop.test[,draw:=draw]
#     pop.test[,model:=model]
#     return(pop.test)
# }
# 
# y<-rbindlist(y)
# 
# save(y, file="/scratch/pae3g/revisions/evolution/single_population_sign_universal_threshold_dropmissing.Rdata")
# 
# 
# 
# 

load("/scratch/pae3g/revisions/evolution/bergland2014_sign_universal_threshold_dropmissing_bottom_quantiles.Rdata")

b2014<-copy(y)

b2014[,permuted:=ifelse(perm==0, F, T)]
b2014[perm!=0&pop=="both", group:="both-permuted"]
b2014[perm==0&pop=="both", group:="both-observed"]
b2014[perm==0&pop=="A", group:="A-observed"]
b2014[perm==0&pop=="B", group:="B-observed"]
b2014[perm!=0&pop=="A", group:="A-permuted"]
b2014[perm!=0&pop=="B", group:="B-permuted"]

b2014[,pheno2:=ifelse(pheno=="diapause.bin", "st. 8", "st. 10")]
b2014[,pheno2:=factor(b2014$pheno2, levels=c("st. 8", "st. 10"))]
b2014[top==-5, th:="Bottom 99.999%"]
b2014[top==-4, th:="Bottom 99.99%"]
b2014[top==-3, th:="Bottom 99.9%"]
b2014[top==-2, th:="Bottom 99%"]
b2014[top==-1, th:="Bottom 90%"]
b2014[top==0, th:="all SNPs"]
b2014[top=="lasso", th:="LASSO"]

b2014[,th:=factor(th, levels=c("all SNPs", "Bottom 90%", "Bottom 99%", "Bottom 99.9%", "Bottom 99.99%", "Bottom 99.999%"))]

b2014.sum<-b2014[,.(med=median(poly), q.025=quantile(poly, 0.025), q.975=quantile(poly, .975)), .(th,group, pheno, pheno2, pop, test,top, permuted)]
b2014.sum[,pheno2:=factor(b2014.sum$pheno2, levels=c("st. 8", "st. 10"))]



stats<-merge(b2014[permuted==F], b2014.sum[permuted==T], by=c("th", "pheno", "pheno2", "pop", "test", "top"))
stats[,over:=poly>q.975]
stats[,under:=poly<q.025]

stats.sum<-stats[,.(n=.N, prop.over=sum(over)/.N*100, prop.under=sum(under)/.N*100, max=max(poly)), .(th, pheno, pheno2, pop, test, top)]

a.plot<-ggplot(b2014.sum[test=="clinal" &top!=(-5)])+
    geom_point(data=b2014.sum[test=="clinal"&top!=(-5)], aes(x=pop, y=med, color=group), position=position_dodge(width=0.5))+
    geom_errorbar(data=b2014.sum[test=="clinal"&top!=(-5)], aes(x=pop, ymax=q.975, ymin=q.025, color=group), width=0.2, position=position_dodge(width=0.5))+
    labs(x="", y="sum(GWAS coefficient*model coefficient)", color="", title="clinal")+
    scale_color_manual(values=c("dodgerblue2", "grey80", "darkorchid", "grey80", "lightseagreen", "grey80"))+
    facet_grid(th~pheno2, scales ="free_y")+
    theme(legend.position = "none")+
    geom_text(data=stats.sum[test=="clinal" & prop.over>50&top!=(-5)], aes(x=pop, y=max+0.1*max, label=prop.over))


# b.plot<-ggplot(b2014.sum[test=="seasonal" ])+
#     geom_point(data=b2014.sum[test=="seasonal"], aes(x=pop, y=med, color=group), position=position_dodge(width=0.5))+
#     geom_errorbar(data=b2014.sum[test=="seasonal"], aes(x=pop, ymax=q.975, ymin=q.025, color=group), width=0.2, position=position_dodge(width=0.5))+
#     labs(x="", y="", color="", title="seasonal")+
#     scale_color_manual(values=c("dodgerblue2", "grey80", "darkorchid", "grey80", "lightseagreen", "grey80"))+
#     facet_grid(th~pheno2, scales ="free_y")+
#     theme(legend.position = "none")+
#     geom_text(data=stats.sum[test=="seasonal" & prop.over>50], aes(x=pop, y=max+0.1*max, label=prop.over))


pdf("/scratch/pae3g/revisions/figures/bergland2014_main_dropmissing_bottomquantiles.pdf", height=8, width=5)
a.plot
dev.off()