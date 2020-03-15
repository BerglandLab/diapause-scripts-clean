
library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(foreach)
library(doMC)
registerDoMC(20)
library(Rmisc)
library(ggbeeswarm)

files<-fread("/scratch/pae3g/genome-reconstruction/universal_input.txt")

ad<-fread("/scratch/pae3g/revisions/evolution/admixture_tracts_ZI.txt")

setnames(ad, c("V1", "V2", "V3", "V4", "V5"), c("line", "CHR", "start", "stop", "length"))
ad[,chr:=tstrsplit(CHR, split="Chr")[[2]]]

admix<-foreach(pop=files$V1[1:3000], phenotype=files$V2[1:3000],draw=files$V3[1:3000], perm=files$V4[1:3000], model=files$V6[1:3000], .errorhandling="remove")%dopar%{
    print(paste(perm, draw, sep=", "))
    load(paste0("/scratch/pae3g/revisions/lasso/", "lasso_", phenotype, "_draw",draw, "_perm", perm, "_", model, "_pop", pop, ".lassoSites.Rdata"))
    sites[, chr:=tstrsplit(site,split="_")[[1]]]
    sites[, pos:=as.numeric(tstrsplit(site,split="_")[[2]])]
    sites<-sites[grep("SNP", site)]
    z<-foreach(chr.test=sites$chr, pos.test=sites$pos)%dopar%{
        t<-ad[chr==chr.test & pos.test>start & pos.test<stop]
        return(data.table(chr=chr.test,
                          pos=pos.test,
                          perm=perm,
                          draw=draw,
                          pheno=phenotype,
                          model=model,
                          pop=pop,
                          n.admix=nrow(t)))
    }
    return(rbindlist(z))
    
}
admix<-rbindlist(admix)

save(admix, file="/scratch/pae3g/revisions/evolution/admix_universal_lasso.Rdata")