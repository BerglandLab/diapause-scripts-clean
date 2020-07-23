
# library(foreach)
# library(data.table)
# library(ggplot2)
# library(cowplot)
# theme_set(theme_cowplot())
# library(doMC)
# registerDoMC(10)
# library(ggbeeswarm)




#try making averaged qqplot
#largemem 16 core


library(foreach)
library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(doMC)
registerDoMC(16)
library(ggbeeswarm)



files<-fread("/scratch/pae3g/genome-reconstruction/universal_input.txt") #on dryad
setnames(files, c("pop", 'pheno', "draw", "perm", "seed", "model"))
files[,permuted:=ifelse(perm==0, F, T)]
setkey(files, pop, pheno, model, permuted)

x<-unique(files[,.(pop, pheno, model, permuted)])

#files[,index:=c(1:nrow(files))]
q<-foreach(i=c(1:nrow(x)))%do%{
    files2=files[J(x[i])]
    #first load in one file to choose the random thinning samples
    load(paste("/scratch/pae3g/revisions/genesis_", files2$pheno[1], "_draw", files2$draw[1], "_perm", files2$perm[1], "_pop", files2$pop[1], "_" , files2$model[1], "_allsnpgrm_wolbachia_dropmissing.Rdat", sep=""))
    a<-copy(assoc.results)
    #a<-a[!is.na(Score.pval)]
    
    #assoc.results[,maf:=pmin(freq, 1-freq)]
   # assoc.results<-assoc.results[maf>=0.05]
    a[,e := -log10( 1:length(Score.pval)/length(Score.pval))]
    #assoc.results[e<2, e.round:=round(e,digits=7)]
       # assoc.results[e>=2, e.round:=e]

    thinned<-sample(a[e<2, e], 10000)
    samples=c(thinned, a[e>=2,e])
    
    p<-foreach(pop=files2$pop, phenotype=files2$pheno, draw=files2$draw, perm=files2$perm, model=files2$model, .errorhandling="remove")%dopar%{
        load(paste("/scratch/pae3g/revisions/genesis_", phenotype, "_draw", draw, "_perm", perm, "_pop", pop, "_" , model, "_allsnpgrm_wolbachia_dropmissing.Rdat", sep=""))
       # print(index)
        #assoc.results[,maf:=pmin(freq, 1-freq)]
        #assoc.results<-assoc.results[maf>=0.05]
        assoc.results[,pop:=pop]
        assoc.results[,draw:=draw]
        assoc.results[,e := -log10( 1:length(Score.pval)/length(Score.pval))]
        assoc.results<-assoc.results[!is.na(Score.pval)]
        assoc.results[,o := -log10(sort(Score.pval,decreasing=F))]
        
        #assoc.results[e>=2, e.round:=e]
        #assoc.results[e<2, e.round:=round(e,7)]
        setkey(assoc.results, e)
        return(assoc.results[J(samples), .(o,e, phenotype, draw, perm, pop, model)])
    }
    p<-rbindlist(p)
    #p[e<1,e.round:=round(e, digits=6)]
    #p[e>=1,e.round:=round(e, digits=2)]
    
    p[, permuted:=ifelse(perm==0, F, T)]
    qq.sum<-p[,.(avg.o=mean(o, na.rm=T), sd.o=sd(o, na.rm=T), n=.N), .(e, phenotype, pop, model, permuted )]
    qq.sum[,se:=sd.o/sqrt(n)]
    return(qq.sum)
}

q<-rbindlist(q)
save(q, file="/scratch/pae3g/revisions/avg_qqplot_data_thinned.Rdat")

load("/scratch/pae3g/revisions/avg_qqplot_data_thinned.Rdat")
q<-q[!is.na(phenotype)&!is.na(model)]

q[permuted==T&pop=="both", group:="both-permuted"]
q[permuted==F&pop=="both", group:="both-observed"]
q[permuted==T&pop=="A", group:="A-permuted"]
q[permuted==F&pop=="A", group:="A-observed"]
q[permuted==T&pop=="B", group:="B-permuted"]
q[permuted==F&pop=="B", group:="B-observed"]

q[phenotype=="diapause.bin9", phenotype:="Stage 10"]
q[phenotype=="diapause.bin", phenotype:="Stage 8"]
q[,phenotype:=factor(phenotype, levels=c("Stage 8", "Stage 10"))]
q[,model:=ifelse(model=="nonloco", "non-LOCO" , "LOCO")]

pdf("/scratch/pae3g/revisions/figures/loco_nonloco_qqplot.pdf", height=8, width=8)
ggplot(q)+
    geom_line(data=q, aes(x=e, y=avg.o, color=group))+
    geom_ribbon(data=q, aes(x=e, ymin=avg.o-sd.o, ymax=avg.o+sd.o, fill=group, alpha=0.3))+
    scale_color_manual(values=c("dodgerblue2", "grey80", "darkorchid",  "grey80", "lightseagreen","grey80")) +
    scale_fill_manual(values=c("dodgerblue2", "grey80", "darkorchid","grey80",  "lightseagreen","grey80")) +
    labs(x=expression("-log"[10]*"(expected P)"),y=expression("-log"[10]*"(observed P)"))+
    facet_grid(phenotype~model)+
    theme(legend.position = "none")+
    geom_abline(slope=1, intercept=0)
    
dev.off()



