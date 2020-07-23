
library(data.table)
library(foreach)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

library(doMC)
registerDoMC(20)

#read in a file that has ZI allele freqs

b<-fread("/scratch/pae3g/revisions/evolution/fst_ZI_AT_gr_12_fall.txt")
#loop through permutations and test pro-diapause alleles at different thresholds

files<-fread("/scratch/pae3g/genome-reconstruction/universal_input.txt")


zambia<-foreach(pop=files$V1[1:3000], phenotype=files$V2[1:3000],draw=files$V3[1:3000], perm=files$V4[1:3000], model=files$V6[1:3000], .errorhandling="remove")%dopar%{
    print(paste(draw, perm, sep=","))
    #read gwas adn do some fixes
    load(paste("/scratch/pae3g/revisions/genesis_", phenotype, "_draw", draw, "_perm", perm, "_pop", pop, "_" , model, "_allsnpgrm_wolbachia_dropmissing.Rdat", sep=""))
    gwas<-assoc.results   
    gwas[,maf:=pmin(freq, 1-freq)]
    gwas<-gwas[maf>=0.05]
    gwas[,q:=log10(frank(Score.pval)/(length(Score.pval)+1))]
    g<-merge(gwas, b, by=c("chr", "pos"))
    g[,pro.diapause.zambia:=ifelse(sign(Score)==1, p2, 1-p2)]
    
    a<-foreach (top = seq(from=-5, to=0, by=1)) %do% {
        
        
        g.sum<-g[q<=top,.(med=median(pro.diapause.zambia, na.rm=T),
                          n=.N,
                          prop.01=sum(pro.diapause.zambia>.01, na.rm=T)/.N,
                          prop.05=sum(pro.diapause.zambia>.05, na.rm=T)/.N,
                          prop.1=sum(pro.diapause.zambia>.1, na.rm=T)/.N,
                          prop.2=sum(pro.diapause.zambia>.2, na.rm=T)/.N)]
        g.sum[,perm:=perm]
        g.sum[,draw:=draw]
        g.sum[,top:=top]
        g.sum[,pop:=pop]
        g.sum[,pheno:=phenotype]
        g.sum[,model:=model]
        return(g.sum)
    }
    return(rbindlist(a))
}


zambia<-rbindlist(zambia)

save(zambia, file="/scratch/pae3g/revisions/evolution/zambia_universal_threshold_dropmissing.Rdata")


# library(data.table)
# library(foreach)
# library(ggplot2)
# library(cowplot)
# theme_set(theme_cowplot())
# library(ggbeeswarm)
# load("/scratch/pae3g/revisions/evolution/zambia_universal_threshold.Rdata")
# 
# lasso.zi.sum.melt<-melt(zambia, id.vars=c("pop", "pheno", "draw", "perm", "top"), measure.vars=c("med", "prop.01", "prop.05", "prop.1", "prop.2"))
# 
# lasso.zi.sum.melt[perm!=0&pop=="both", group:="both-permuted"]
# lasso.zi.sum.melt[perm==0&pop=="both", group:="both-observed"]
# lasso.zi.sum.melt[perm==0&pop=="A", group:="A-observed"]
# lasso.zi.sum.melt[perm==0&pop=="B", group:="B-observed"]
# lasso.zi.sum.melt[perm!=0&pop=="A", group:="A-permuted"]
# lasso.zi.sum.melt[perm!=0&pop=="B", group:="B-permuted"]
# lasso.zi.sum.melt[,pheno2:=ifelse(pheno=="diapause.bin", "stage 8", "stage 10")]
# 
# lasso.zi.sum.melt[top==-5, th:="Top 0.001%"]
# lasso.zi.sum.melt[top==-4, th:="Top 0.01%"]
# lasso.zi.sum.melt[top==-3, th:="Top 0.1%"]
# lasso.zi.sum.melt[top==-2, th:="Top 1%"]
# lasso.zi.sum.melt[top==-1, th:="Top 10%"]
# lasso.zi.sum.melt[top==0, th:="all SNPs"]
# 
# y.sum<-lasso.zi.sum.melt[,.(med=median(value), q.025=quantile(value, 0.025), q.975=quantile(value, .975)), .(group, pheno, pheno2, th,  pop,variable)]
# 
# 
# y.sum[,facet:=paste(pheno2, th, sep=": ")]
# y.sum[,facet:=factor(y.sum$facet, levels=unique(y.sum$facet)[c(1,7,2,8,3,9,4,10,5,11,6,12)])]
# 
# a.plot<-ggplot(y.sum)+
#     geom_point(data=y.sum, aes(x=variable, y=med, color=group), position=position_dodge(width=0.75))+
#     geom_errorbar(data=y.sum, aes(x=variable, ymax=q.975, ymin=q.025, color=group), width=0.2, position=position_dodge(width=0.75))+
#     labs(x="", y="Frequency or Proportion",  color="")+
#     scale_x_discrete(labels=c("median freq", "prop > 0.01", "prop > 0.05", "prop > 0.1", "prop > 0.2"))+
#     scale_color_manual(values=c("#39568CFF", "grey20", "#440154FF", "grey50", "lightseagreen", "grey80"))+
#     theme(legend.position = "none")+
#     facet_wrap(~facet, nrow=6, scales="free_y")+
#     theme(axis.text.x=element_text(angle=45,hjust=1))
# 
# 
# pdf("/scratch/pae3g/revisions/figures/zambia_threshold.pdf", height=8, width=8)
# a.plot
# dev.off()
# 
# 
# 
# 
# 
# zi.plot<-ggplot(lasso.zi.sum.melt[top!=(-5)&top!=0], aes(x=variable, y=value, color=group))+
#     geom_beeswarm(dodge.width=.8,alpha=0.25, size=0.05)+
#     labs(x="", y="freq. or proportion", color="")+
#     scale_color_manual(values=c("#39568CFF", "#440154FF", "lightseagreen", "grey80"))+
#     scale_x_discrete(labels=c("median", "> 0.01", "> 0.05", "> 0.10", "> 0.20"))+
#     theme(strip.background = element_blank())+
#     facet_grid(pheno2~th,scales="free")+
#     theme(legend.position = "none")+
#     theme(axis.text.x=element_text(angle=45, hjust=1))
# 
# pdf("/scratch/pae3g/revisions/figures/zambia_threshold.pdf", height=8, width=10)
# zi.plot
# dev.off()
# 
# #calculate statistics
# 
# foreach(th=c(-5:0))%do%{
#     t<-zambia[top==th]
#     
# }
# 
