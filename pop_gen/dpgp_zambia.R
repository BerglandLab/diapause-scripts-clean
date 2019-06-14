
#are more pro-diapause alleles segregating in AFrica than we would expect by chance?

# take different thresholds, calculate the fraction of pro-diapause alleles > 0.05 frequency, and compare to permutations

library(data.table)
library(foreach)
library(cowplot)
library(ggbeeswarm)
library(doMC)
registerDoMC(3)

#read in a file that has ZI allele freqs

    b<-fread("/mnt/pricey_2/priscilla/dpgp/fst_ZI_AT_gr_12_fall.txt")
    #loop through permutations and test pro-diapause alleles at different thresholds
    
    perm.test<-foreach(perm=c(0, 101:200))%do%{
        print(perm)
        gwas<-fread(paste0("/mnt/sammas_storage/bergland-lab/Priscilla/perm_summaries/gwas_p_score_inv_id_perm", perm, ".txt"))
        gwas<-merge(gwas, b, by=c("chr", "pos"))
        gwas[,pro.diapause:=ifelse(sign(gwas.Score)==1, p2, 1-p2)]
        
        g<-foreach(top=c(101, 490, 5316, 25250, 67914))%dopar%{
            snps<-fread(paste0("/mnt/sammas_storage/bergland-lab/Priscilla/perm_summaries/perm", perm, "_top", top, "_ld0.4_200kb.clumped"))
            setnames(snps, c("CHR", "BP"), c("chr", "pos"))
            snps[,chr:=as.character(chr)]
            snps[chr=="23", chr:="X"]
            snps=merge(snps[,.(chr, pos)], gwas, by=c("chr", "pos"))
            snps[,perm:=perm]
            snps[,top:=top]
            return(snps)
        }
        return(rbindlist(g))
    }
    
    perm.test<-rbindlist(perm.test)
    
    write.table(perm.test, "/mnt/pricey_2/priscilla/dpgp/pro_diapause_alleles_africa_clump200kb.txt", quote=F, sep="\t", row.names = F)
    
    perm.test<-fread("/mnt/pricey_2/priscilla/dpgp/pro_diapause_alleles_africa_clump200kb.txt")
    
    #calculate the proportion of SNPs present above given allele frequencies
    perm.test.sum<-perm.test[,.(n=.N,
                                prop.01=sum(pro.diapause>0.01, na.rm=T)/.N,
                                prop.05=sum(pro.diapause>0.05, na.rm=T)/.N,
                                prop.1=sum(pro.diapause>0.1, na.rm=T)/.N,
                                prop.2=sum(pro.diapause>0.2, na.rm=T)/.N,
                                med.p=median(pro.diapause, na.rm=T),
                                max.p=max(pro.diapause, na.rm=T),
                                min.p=min(pro.diapause, na.rm=T)),
                             .(top, perm)]
    #melt this data table
    perm.test.melt<-melt(perm.test.sum, id.vars=c("top", "perm"), measure.vars=c("prop.01", "prop.05", "prop.1", "prop.2", "med.p"))

    #calculate .05 quantile of permutations for one-tailed test (prediction is deficit of SNPs)
    perm.test.q<-perm.test.sum[perm!=0, .(q.05.prop.01=quantile(prop.01, .05),
                                          q.05.prop.05=quantile(prop.05, .05),
                                          q.05.prop.1=quantile(prop.1, .05),
                                          q.05.prop.2=quantile(prop.2, .05),
                                          q.05.med.p=quantile(med.p, .05)), 
                               .(top)]    
    #merge quantiles with actual data
    perm.test.q<-merge(perm.test.sum[perm==0], perm.test.q, by="top")
    
    #see which ones pass
    perm.test.q[prop.01<=q.05.prop.01]
    perm.test.q[prop.05<=q.05.prop.05]
    perm.test.q[prop.1<=q.05.prop.1]
    perm.test.q[prop.2<=q.05.prop.2]
    perm.test.q[med.p<=q.05.med.p]
    
