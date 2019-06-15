#get top snps from permutations for haplotype block calculation

library(data.table)
library(foreach)
a<-fread("/mnt/sammas_storage/bergland-lab/Priscilla/perm_summaries/gwas_p_score_inv_id_perm0.txt")
a[,fdr:=p.adjust(gwas.p, method="fdr")]

#count the number of snps below fdr cutofs
n1=length(a[fdr<.025, id])
n2=length(a[fdr<.05, id])
n3=length(a[fdr<.1, id])
n4=length(a[fdr<.15, id])
n5=length(a[fdr<.2, id])


#go through permutations to calculate the p-value that corresponds to the same # of snps as in the FDR cutoff
pt<-foreach(perm=c(0,101:200))%do%{
    print(perm)
    p<-fread(paste0("/mnt/sammas_storage/bergland-lab/Priscilla/perm_summaries/gwas_p_score_inv_id_perm", perm, ".txt"))
    #make a snp id column for plink
    p[,SNP:=paste(chr, pos, "SNP", sep="_")]
    #save a table to be read by plink
    write.table(p[,.(SNP, gwas.p)], paste0("/mnt/sammas_storage/bergland-lab/Priscilla/perm_summaries/perm", perm, "_gwas_for_clump.txt"), quote=F, sep="\t", row.names=F)
    #order by ascending pvalue
    p<-p[order(gwas.p)]
    return(data.table(perm=perm,
                      fdr=c(0.025, 0.05, 0.1, 0.15, 0.2),
                      top=c(n1, n2, n3, n4, n5),
                      #pull out pvalues at each row number
                      pthresh=format(c(p[n1, gwas.p], p[n2, gwas.p], p[n3, gwas.p], p[n4, gwas.p], p[n5, gwas.p]), scientific=F)))
    }
pt<-rbindlist(pt)

#save parameters to use in bash script
write.table(pt, "/mnt/pricey_2/priscilla/clump_p_thresholds_perms.txt", row.names=F, col.names=F, quote=F, sep="\t")




