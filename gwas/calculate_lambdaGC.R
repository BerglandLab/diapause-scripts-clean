library(data.table)
library(foreach)
library(doMC)
registerDoMC(16)
library(cowplot)
library(ggbeeswarm)

y<-foreach(perm=c(0, 101:200), .combine="rbind")%dopar%{
    print(perm)
    a<-fread(paste0("/scratch/pae3g/evolution/gwas_p_score_inv_id_perm", perm, ".txt"))
    l<-median(qchisq(1-a$gwas.p,1))/qchisq(0.5,1)
    l.2L<-median(qchisq(1-a[chr=="2L"]$gwas.p,1))/qchisq(0.5,1)
    l.2R<-median(qchisq(1-a[chr=="2R"]$gwas.p,1))/qchisq(0.5,1)
    l.3L<-median(qchisq(1-a[chr=="3L"]$gwas.p,1))/qchisq(0.5,1)
    l.3R<-median(qchisq(1-a[chr=="3R"]$gwas.p,1))/qchisq(0.5,1)
    l.X<-median(qchisq(1-a[chr=="X"]$gwas.p,1))/qchisq(0.5,1)
    
    return(data.table(perm=perm,
                      lambda=l,
                      l.2L=l.2L,
                      l.2R=l.2R,
                      l.3L=l.3L,
                      l.3R=l.3R,
                      l.X=l.X))
}
y<-rbindlist(y)
write.table(y, "/scratch/pae3g/evolution/lambda.txt", quote=F, sep="\t", row.names = F)