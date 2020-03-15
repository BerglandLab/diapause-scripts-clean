#look at annotations of top SNPs

#first look at lasso

library(data.table)
library(foreach)


load("/scratch/pae3g/revisions/snp_annotations.Rdat")

load("/scratch/pae3g/revisions/all_lasso_for_alan.Rdata")
m[col.class%in%c("3_prime_UTR_variant", "5_prime_UTR_premature_start_codon_gain_variant", "5_prime_UTR_variant"), group:="UTR"]
m[col.class%in%c("downstream_gene_variant", "upstream_gene_variant"), group:="up/down"]
m[col.class=="synonymous_variant"|col.class=="splice_region_variant&synonymous_variant", group:="syn"]
m[col.class=="missense_variant"|col.class=="missense_variant&splice_region_variant", group:="non-syn"]
m[col.class=="intergenic_region", group:="intergenic"]
m[col.class%in%c("splice_region_variant&intron_variant","intron_variant", "splice_acceptor_variant&intron_variant", "splice_donor_variant&intron_variant"), group:="intron"]

p<-merge(p, m, by=c("chr", "pos"))

p.total<-p[,.(n.total=.N), .(perm, GRM, pop, pheno)]

class.test<-foreach(type=unique(m$group))%do%{
    a<-p[group==type]
    a.sum<-a[,.(n=.N), .(perm, GRM, pop, pheno)]
    a.sum[,group:=type]
    return(a.sum)
}

class.test<-rbindlist(class.test)


test<-foreach(test.type=unique(m$group))%do%{
    b<-merge(p.total, class.test[group==test.type], by=c("GRM", "perm", "pheno", "pop"), all.x=T)
    b[is.na(n), n:=0]
    b[,prop:=n/n.total]
    x<-foreach(pop.test=c("A", "B", "both"))%do%{
        y<-foreach(pheno.test=c("diapause.bin", "diapause.bin9"))%do%{
            med.obs<-median(b[pheno==pheno.test&pop==pop.test&perm==0, prop])
            func<-ecdf(b[pheno==pheno.test&pop==pop.test&perm!=0,prop])
            prob<-func(med.obs)
            return(data.table(test=test.type,
                              pop=pop.test,
                              pheno=pheno.test,
                              p=prob,
                              q="LASSO",
                              med=med.obs))
        }
        return(rbindlist(y))
    }
    return(rbindlist(x))
}
        
test<-rbindlist(test)

test<-test[!is.na(test)]

#look at top percentages
load("/scratch/pae3g/revisions/gwas_top1percent.Rdat")
setnames(y, "draw", "GRM")
setnames(y, "phenotype", "pheno")
#caculate number of snps in each quantile
q.total<-foreach(q.t=c(-4, -3, -2))%do%{
    z<-y[q<=q.t,.(n.total=.N), .(perm, GRM, pop, pheno)]
    z[,q:=q.t]
    return(z)
}
q.total<-rbindlist(q.total)

#merge wtih annotations
y<-merge(y, m, by=c("chr", "pos"))

#count number of each type in each gwas
class.test.q<-foreach(type=unique(m$group))%do%{
    q.b<-foreach(q.t=c(-4, -3, -2))%do%{
        a<-y[group==type&q<=q.t]
        a.sum<-a[,.(n=.N), .(perm, GRM, pop, pheno)]
        a.sum[,group:=type]
        a.sum[,q:=q.t]
        return(a.sum)
        }
    return(rbindlist(q.b))
}

class.test.q<-rbindlist(class.test.q)

test.q<-foreach(q.t=c(-4, -3, -2))%do%{
    test.z<-foreach(test.type=unique(m$group))%do%{
        b<-merge(q.total[q==q.t], class.test.q[group==test.type&q==q.t], by=c("GRM", "perm", "pheno", "pop"), all.x=T)
        b[is.na(n), n:=0]
        b[,prop:=n/n.total]
        x<-foreach(pop.test=c("A", "B", "both"))%do%{
            y<-foreach(pheno.test=c("diapause.bin", "diapause.bin9"))%do%{
                med.obs<-median(b[pheno==pheno.test&pop==pop.test&perm==0, prop])
                func<-ecdf(b[pheno==pheno.test&pop==pop.test&perm!=0,prop])
                prob<-func(med.obs)
                return(data.table(test=test.type,
                                  pop=pop.test,
                                  pheno=pheno.test,
                                  p=prob,
                                  q=q.t,
                                  med=med.obs))
            }
            return(rbindlist(y))
        }
        return(rbindlist(x))
    }
    
    test.z<-rbindlist(test.z)
    return(test.z)
}

test.q<-rbindlist(test.q)

test.q<-test.q[!is.na(test)]

annotation.test<-rbind(test, test.q)
annotation.test.wide<-dcast(annotation.test,test+pop+pheno~q , value.var="p")
annotation.test.wide.med<-dcast(annotation.test,test+pop+pheno~q , value.var="med")

annotation.test.wide[,phenotype:=ifelse(pheno=="diapause.bin", 'stage 8', "stage 10")]
annotation.test.wide[,pheno:=NULL]
setnames(annotation.test.wide, c("-2", "-3", "-4"), c("Top 1%", "Top 0.1%", "Top 0.01%"))
write.csv(annotation.test.wide, "/scratch/pae3g/revisions/figures/annotation_enrichment.csv")



annotation.test.wide.med[,phenotype:=ifelse(pheno=="diapause.bin", 'stage 8', "stage 10")]
annotation.test.wide.med[,pheno:=NULL]
setnames(annotation.test.wide.med, c("-2", "-3", "-4"), c("Top 1%", "Top 0.1%", "Top 0.01%"))
write.csv(annotation.test.wide.med, "/scratch/pae3g/revisions/figures/annotation_proportions.csv")