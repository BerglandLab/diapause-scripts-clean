
#are more pro-diapause alleles segregating in AFrica than we would expect by chance?

# take different thresholds, calculate the fraction of pro-diapause alleles > 0.05 frequency, and compare to permutations

library(data.table)
library(foreach)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(ggbeeswarm)
library(doMC)
registerDoMC(20)

files<-fread("/scratch/pae3g/genome-reconstruction/universal_input.txt")

ad<-fread("/scratch/pae3g/revisions/evolution/admixture_tracts_ZI.txt")

setnames(ad, c("V1", "V2", "V3", "V4", "V5"), c("line", "CHR", "start", "stop", "length"))
ad[,chr:=tstrsplit(CHR, split="Chr")[[2]]]
ad[,id:=1:nrow(ad)]
admix<-foreach(pop=files$V1[1:3000], phenotype=files$V2[1:3000],draw=files$V3[1:3000], perm=files$V4[1:3000], model=files$V6[1:3000])%dopar%{
    print(paste(draw, perm, sep=","))
    #read lasso SNPs and do some fixes
    load(paste("/scratch/pae3g/revisions/genesis_", phenotype, "_draw", draw, "_perm", perm, "_pop", pop, "_" , model, "_allsnpgrm_wolbachia.Rdat", sep=""))
    gwas<-assoc.results
    gwas[,maf:=pmin(freq, 1-freq)]
    gwas<-gwas[maf>=0.05]
    gwas[,q:=log10(frank(Score.pval)/(length(Score.pval)+1))]
    setkey(gwas, chr, pos)
    z<-foreach(chr.test=ad$chr, start.test=ad$start, stop.test=ad$stop, id.test=ad$id)%do%{
        t<-gwas[chr==chr.test & pos>start.test & pos<stop.test]
        t[,id:=id.test]
        return(t)
    }
    z<-rbindlist(z)
    a<-foreach (top = seq(from=-5, to=0, by=1)) %do% { 
        w.sum<-z[q<=top, .(total.admix=.N, unique.admix=length(unique(variant.id)))]
        w.sum[,perm:=perm]
        w.sum[,draw:=draw]
        w.sum[,top:=top]
        w.sum[,pheno:=phenotype]
        w.sum[,pop:=pop]
        w.sum[,model:=model]
        w.sum[,ntotal:=nrow(gwas[q<=top])]
        return(w.sum)

    }
    return(rbindlist(a))
}

admix<-rbindlist(admix)

save(admix, file="/scratch/pae3g/revisions/evolution/ZI_admix_universal_threshold.Rdata")

# admix.sum<-admix[,.(total.admix=sum(n.admix)/.N, unique.admix=sum(n.admix>0)/.N, n=.N), .(perm, draw)]
# 
# 
# a<-ggplot(admix.sum, aes(x=total.admix, y=..scaled.., fill=ifelse(perm==0, "obs", "perm")))+
#     geom_density(alpha=0.5)+
#     scale_fill_manual(values=c("lightseagreen", "grey80"))+
#     labs(x="total # admixture overlaps / total # SNPs ")+
#     theme(legend.position = "none")
# 
# b<-ggplot(admix.sum, aes(x=unique.admix, y=..scaled.., fill=ifelse(perm==0, "obs", "perm")))+
#     geom_density(alpha=0.5)+
#     scale_fill_manual(values=c("lightseagreen", "grey80"))+
#     labs(x="# of unique SNPs in admixture regions / total # of SNPs")+
#     theme(legend.position = "none")
# 
# plot_grid(a,b, nrow=1 )