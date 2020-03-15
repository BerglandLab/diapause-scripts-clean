
library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(foreach)
library(doMC)
registerDoMC(16)

load("/scratch/pae3g/revisions/all_lasso_for_alan.Rdata")
p<-p[chr%in%c("2L", "2R", "3L", "3R", "X")]
p[,GRM:=as.integer(GRM)]
p[,perm:=as.integer(perm)]


#get list of all samples
files<-fread("/scratch/pae3g/genome-reconstruction/universal_input.txt")
setnames(files, c("pop", "pheno", "GRM", "perm", "seed", "model"))
files<-files[model=="nonloco"]
files[,GRM:=as.integer(GRM)]

#A and B overlaps in LASSO
#get totals
p.n<-p[,.(n=.N), .(perm, GRM, pheno, pop)]

#overlap lasso snps in A and B, keeping track of permutation
a.b<-merge(p[pop=="A"], p[pop=="B"], by=c("chr", 'pos', "perm", "GRM", "pheno"))
a.b.n<-a.b[,.(n=.N), .(perm, GRM, pheno)]
a.b.m<-merge(files[pop=="A",.(GRM, perm, pheno)], a.b.n,  all.x=T, by=c("pheno", "GRM", "perm"))
a.b.m[,q:="lasso"]
a.b.m[,type:="a-b"]

a.both<-merge(p[pop=="A"], p[pop=="both"], by=c("chr", 'pos', "perm", "GRM", "pheno"))
a.both.n<-a.both[,.(n=.N), .(perm, GRM, pheno)]
a.both.m<-merge(files[pop=="A",.(GRM, perm, pheno)], a.both.n,  all.x=T, by=c("pheno", "GRM", "perm"))
a.both.m[,q:="lasso"]
a.both.m[,type:="a-both"]


b.both<-merge(p[pop=="B"], p[pop=="both"], by=c("chr", 'pos', "perm", "GRM", "pheno"))
b.both.n<-b.both[,.(n=.N), .(perm, GRM, pheno)]
b.both.m<-merge(files[pop=="A",.(GRM, perm, pheno)], b.both.n,  all.x=T, by=c("pheno", "GRM", "perm"))
b.both.m[,q:="lasso"]
b.both.m[,type:="b-both"]

lasso.pop<-rbind(a.b.m, a.both.m, b.both.m)
#stage 7 and stage 9 overlaps in lasso
d7.d9<-merge(p[pheno=="diapause.bin9"], p[pheno=="diapause.bin"], by=c("chr", 'pos', "perm", "GRM", "pop"))
d7.d9.n<-d7.d9[,.(n=.N), .(perm, GRM, pop)]

pheno.lasso<-merge(d7.d9.n, files[pheno=="diapause.bin9", .(GRM, perm, pop)], all.y=T, by=c("GRM", "perm", "pop"))
pheno.lasso[,q:="lasso"]


#now read in GWAS, take top 1% with quantile rankings and return
# files<-fread("/scratch/pae3g/genome-reconstruction/universal_input.txt")
# 
# y<-foreach(pop=files$V1[1:3000], phenotype=files$V2[1:3000],draw=files$V3[1:3000], perm=files$V4[1:3000], model=files$V6[1:3000], .errorhandling="remove")%dopar%{
#     print(paste(draw, perm, sep=","))
#     #read gwas adn do some fixes
#     load(paste("/scratch/pae3g/revisions/genesis_", phenotype, "_draw", draw, "_perm", perm, "_pop", pop, "_" , model, "_allsnpgrm_wolbachia.Rdat", sep=""))
#     gwas<-assoc.results
#     gwas[,maf:=pmin(freq, 1-freq)]
#     gwas<-gwas[maf>=0.05]
#     gwas[,q:=log10(frank(Score.pval)/(length(Score.pval)+1))]
#     gwas[,pop:=pop]
#     gwas[,draw:=draw]
#     return(gwas[q<=(-2)])
# }
# 
# y<-rbindlist(y)
# save(y, file="/scratch/pae3g/revisions/gwas_top1percent.Rdat")

load("/scratch/pae3g/revisions/gwas_top1percent.Rdat")
files<-fread("/scratch/pae3g/genome-reconstruction/universal_input.txt")
setnames(files, c("pop", "phenotype", "GRM", "perm", "seed", "model"))
files<-files[model=="nonloco"]
files[,draw:=as.integer(GRM)]
files[,perm:=as.character(perm)]

setkey(y, q, pop)
pop.overlap<-foreach(q.th=c(-4, -3, -2))%do%{
    g<-y[q<=q.th]
    a.b<-merge(g[pop=="A"], g[pop=="B"], by=c("chr", 'pos', "perm", "draw", "phenotype"))
    a.b.n<-a.b[,.(n=.N), .(perm, draw, phenotype)]
    a.b.m<-merge(files[pop=="A",.(draw, perm, phenotype)], a.b.n,  all.x=T, by=c("phenotype", "draw", "perm"))
    a.b.m[,type:="a-b"]
    a.b.m[,q:=q.th]


    a.both<-merge(g[pop=="A"], g[pop=="both"], by=c("chr", 'pos', "perm", "draw", "phenotype"))
    a.both.n<-a.both[,.(n=.N), .(perm, draw, phenotype)]
    a.both.m<-merge(files[pop=="A",.(draw, perm, phenotype)], a.both.n,  all.x=T, by=c("phenotype", "draw", "perm"))
    a.both.m[,type:="a-both"]
    a.both.m[,q:=q.th]



    b.both<-merge(g[pop=="B"], g[pop=="both"], by=c("chr", 'pos', "perm", "draw", "phenotype"))
    b.both.n<-b.both[,.(n=.N), .(perm, draw, phenotype)]
    b.both.m<-merge(files[pop=="A",.(draw, perm, phenotype)], b.both.n,  all.x=T, by=c("phenotype", "draw", "perm"))
    b.both.m[,type:="b-both"]

    b.both.m[,q:=q.th]

    return(rbind(a.b.m, a.both.m, b.both.m))
}
pop.overlap<-rbindlist(pop.overlap)
setnames(lasso.pop, c("pheno", "GRM"), c("phenotype", "draw"))
pop.overlap<-rbind(pop.overlap, lasso.pop)
pop.overlap[is.na(n), n:=0]
pop.overlap[,group:=ifelse(perm==0, "observed", "permuted")]



#     #stage 7 and stage 9 overlaps in lasso
pheno.overlaps<-foreach(q.th=c(-4, -3, -2))%do%{
    g<-y[q<=q.th]
    
    d7.d9<-merge(g[phenotype=="diapause.bin9"], g[phenotype=="diapause.bin"], by=c("chr", 'pos', "perm", "draw", "pop"))
    d7.d9.n<-d7.d9[,.(n=.N), .(perm, draw, pop)]

    pheno.o<-merge(d7.d9.n, files[phenotype=="diapause.bin9", .(draw, perm, pop)], all.y=T, by=c("draw", "perm", "pop"))
    pheno.o[,q:=q.th]
    return(pheno.o)
    }


pheno.overlaps<-rbindlist(pheno.overlaps)
setnames(pheno.lasso, "GRM", "draw")
pheno.overlaps<-rbind(pheno.overlaps, pheno.lasso)

pheno.overlaps[perm!=0&pop=="both", group:="both-permuted"]
pheno.overlaps[perm==0&pop=="both", group:="both-observed"]
pheno.overlaps[perm==0&pop=="A", group:="A-observed"]
pheno.overlaps[perm==0&pop=="B", group:="B-observed"]
pheno.overlaps[perm!=0&pop=="A", group:="A-permuted"]
pheno.overlaps[perm!=0&pop=="B", group:="B-permuted"]
pheno.overlaps[is.na(n), n:=0]

#summarize
pop.sum<-pop.overlap[,.(med=median(n), q025=quantile(n, .025, na.rm=T), q975=quantile(n, .975)), .(type, phenotype, q, group)]
pop.sum[,pheno2:=ifelse(phenotype=="diapause.bin", "stage 8", "stage 10")]
pop.sum[,pheno2:=factor(pop.sum$pheno2, levels=c("stage 8", "stage 10"))]

pop.sum[q==-4, th:="Top 0.01%"]
pop.sum[q==-3, th:="Top 0.1%"]
pop.sum[q==-2, th:="Top 1%"]

pop.sum[q=="lasso", th:="LASSO"]
pheno.sum<-pheno.overlaps[,.(med=median(n), q025=quantile(n, .025, na.rm=T), q975=quantile(n, .975)), .( pop, q, group)]

pheno.sum[q==-4, th:="Top 0.01%"]
pheno.sum[q==-3, th:="Top 0.1%"]
pheno.sum[q==-2, th:="Top 1%"]

pheno.sum[q=="lasso", th:="LASSO"]

ab.plot<-ggplot(pop.sum[type=="a-b"])+
    geom_point(data=pop.sum[type=="a-b"], aes(x=pheno2, y=med, color=group), position=position_dodge(width=0.5))+
    geom_errorbar(data=pop.sum[type=="a-b"], aes(x=pheno2, ymax=q975, ymin=q025, color=group), width=0.2, position=position_dodge(width=0.5))+
    labs(x="", y="number of overlapping SNPs between populations", color="")+
    scale_color_manual(values=c("palevioletred2", "grey80"))+
    facet_grid(th~., scales ="free_y")+
    theme(legend.position="none")+



c.plot<-ggplot(pheno.sum)+
    geom_point(data=pheno.sum, aes(x=pop, y=med, color=group), position=position_dodge(width=0.5))+
    geom_errorbar(data=pheno.sum, aes(x=pop, ymax=q975, ymin=q025, color=group), width=0.2, position=position_dodge(width=0.5))+
    labs(x="", y="number of overlapping SNPs between phenotypes", color="")+
    scale_color_manual(values=c("dodgerblue2", "grey80", "darkorchid", "grey80", "lightseagreen", "grey80"))+
    facet_grid(th~., scales ="free_y")+
    theme(legend.position="none")



pdf("/scratch/pae3g/revisions/figures/number_of_mapping_overlaps.pdf", height=6, width=8)
plot_grid(ab.plot,c.plot, nrow=1,labels=c("A", "B") , rel_widths = c(.4, .6), align="h", axis="b")
dev.off()