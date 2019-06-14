
library(data.table)
library(foreach)
library(cowplot)
library(doMC)
registerDoMC(20)

seas.th.list=seq(from=-5, to=-.8, by=0.1)

pops=c("AGA_14", "LMA_14", "TKA_14")

pop.test<-foreach(pop=pops)%do%{
    print(pop)
    seas<-fread(paste0("/scratch/pae3g/evolution/fisher_exactJ_", pop,".coef_minp2.txt"))
    setnames(seas, "chrom", "chr")
    
    clump.enrich<-foreach(perm =c(0,101:200)) %dopar% { #cycle through original data and 100 permutations
        print(perm)
        gwas<-fread(paste0("/scratch/pae3g/evolution/gwas_p_score_inv_id_perm", perm, ".txt"))
        a<-foreach (top = c(101, 490, 5316, 25250, 67914)) %do% { #these numbers correspond to the top # of snps in our different FDR cutoffs
            
            snps<-fread(paste0("/nv/vol186/bergland-lab/Priscilla/perm_summaries/perm", perm, "_top", top, "_ld0.4_200kb.clumped"))
            setnames(snps, c("CHR", "BP"), c("chr", "pos"))
            snps[,chr:=as.character(chr)]
            snps[chr=="23", chr:="X"]
            gwas.m<-merge(gwas, seas, by=c("chr", "pos"))
            
            #quantile rank and normalize fisher P values in merged datset
            gwas.m[,q.s:=log10(frank(minp2)/(length(minp2)+1))]
            gwas.m<-merge(gwas.m, snps[, .(chr, pos, SNP)], by=c("chr", "pos"))
            nsnps=nrow(gwas.m)
            #return actual overlapping SNPs
            seasonal<-foreach(seas.q.thresh = seas.th.list, .errorhandling="remove")%dopar%{
                x=gwas.m[q.s<=seas.q.thresh, .(SNP)]
                x[,seas.th:=seas.q.thresh]
                x[,pop:=pop]
                x[,perm:=perm]
                x[,gwas.th:=top]
                x[,n.index.merge:=nsnps]
                return(x)
            }
            return(rbindlist(seasonal))
        }
        return(rbindlist(a))
    }
    
    clump.enrich<-rbindlist(clump.enrich)

    return(clump.enrich)
}

pop.test<-rbindlist(pop.test)

write.table(pop.test, "/scratch/pae3g/evolution/3poptest.txt", quote=F, sep="\t", row.names=F)

#how many times does each snp show up at each threshold combination?


#analayze
library(data.table)

pop.test<-fread("/scratch/pae3g/evolution/3poptest.txt")
p.sum.bysnp<-pop.test[,.(n=.N), .(seas.th, gwas.th, perm, SNP)]
p.sum.bysnp[n==3]

#sanity check
# sum(p.sum.bysnp$n) #ok sum is correct

#for each threshold combination, how many snps overlap in all three populations? this is like the middle number for a triple venn diagram for each threshold combination
p.sum.count<-p.sum.bysnp[n==3,.(n.triple=.N), .(perm, gwas.th, seas.th)]

#how many total snps went into the test at each threshold combination, not just the ones that returned SNPs? use the original output of the enrichment test to get these numbers

pops=c("AGA_14", "LMA_14", "TKA_14")

clump.enrich<-fread("/scratch/pae3g/evolution/index_snps_enrich_by_pop_seas2019.txt")
clump.enrich<-clump.enrich[pop%in%c("AGA_14", "LMA_14", "TKA_14")]
clump.enrich[,n.total:=TT+TF]
p.sum.bypop<-clump.enrich[,. (avg.index.snps=mean(n.total), avg.q.snps=mean(TT)), .(gwas.th, perm, seas.th)]

p.sum.bypop<-clump.enrich[,. (avg.index.snps=mean(n.total), 
                              n.aga=n.total[pop=="AGA_14"], 
                              n.tka=n.total[pop=="TKA_14"], 
                              n.lma=n.total[pop=="LMA_14"], 
                              q.aga=TT[pop=="AGA_14"], 
                              q.tka=TT[pop=="TKA_14"], 
                              q.lma=TT[pop=="LMA_14"], 
                              avg.q.snps=mean(TT)), .(gwas.th, perm, seas.th)]

p.sum.bypop[,prop.aga:=q.aga/n.aga]
p.sum.bypop[,prop.tka:=q.tka/n.tka]
p.sum.bypop[,prop.lma:=q.lma/n.lma]
venn<-merge(p.sum.count, p.sum.bypop, by=c("perm", "gwas.th", "seas.th"), all=T)
venn[,expect:=round(prop.aga*prop.tka*prop.lma*avg.index.snps)]
venn[is.na(n.triple), n.triple:=0]

venn[,chi.stat:=((n.triple-expect)^2)/expect]
venn[,chi.p:=pchisq(chi.stat, 1, lower.tail=F)]


venn[,prop.index:=n.triple/avg.index.snps]
venn[,prop.q:=n.triple/avg.q.snps]


venn.sum<-venn[perm!=0, .(q95.prop.index=quantile(prop.index, .95, na.rm=T), q95.prop.q=quantile(prop.q, .95, na.rm=T)), .(seas.th, gwas.th)]

venn.sum<-merge(venn.sum, venn[perm==0], by=c("seas.th", "gwas.th"))
#chi square test on product of probabilities
#p1*p2*p3*number of snps --> calculate chi square statistic


#look at a couple snps that look interesting