#calculate polygenic scores for clinal and seasonal signal across permutations

##maybe not this, try this instead:
library(data.table)
library(foreach)
library(cowplot)
library(doMC)
registerDoMC(10)


clump.enrich<-foreach(perm =c(0,101:200)) %do% { #cycle through original data and 100 permutations
    print(perm)
    a<-foreach (top = c(101, 490, 5316, 25250, 67914)) %do% { #these numbers correspond to the top # of snps in our different FDR cutoffs
    
        gwas<-fread(paste0("/scratch/pae3g/evolution/gwas_p_score_inv_id_perm", perm, ".txt"))
        snps<-fread(paste0("/nv/vol186/bergland-lab/Priscilla/perm_summaries/perm", perm, "_top", top, "_ld0.4_200kb.clumped"))
        cline<-fread("/scratch/pae3g/evolution/east_coast_cline_V2_clean.txt")
        
        setnames(snps, c("CHR", "BP"), c("chr", "pos"))
        snps[,chr:=as.character(chr)]
        snps[chr=="23", chr:="X"]
        gwas<-merge(gwas, cline, by=c("chr", "pos"))
        #correct sign of cline for reference allele
        gwas[,clinal.beta:=clinal.beta*-1]
        gwas[,clinal.beta.q:=frank(clinal.beta)/(length(clinal.beta)+1)]
        gwas[,clinal.beta.Z.q:=qnorm(clinal.beta.q, 0,1)]
        
        #quantile rank and normalize Scores in merged datset
        gwas[,Score.q:=frank(gwas.Score)/(length(gwas.Score)+1)]
        gwas[,Score.Z.q:=qnorm(Score.q, 0,1)]
        #take product of Score * clinal beta
        gwas[,Score.cline.beta:=gwas.Score*clinal.beta]
        gwas[,Score.Z.q.cline.beta:=Score.Z.q*clinal.beta]
        gwas[,Score.Z.q.cline.beta.Z.q:=Score.Z.q*clinal.beta.Z.q]
        gwas[,Score.cline.beta.Z.q:=gwas.Score*clinal.beta.Z.q]
        
        
        #test for concordance of cline:
        gwas[,conc:=sign(clinal.beta)==sign(gwas.Score)]
        #take product of Score.Z*clinal.beta
        gwas<-merge(gwas, snps[, .(chr, pos)], by=c("chr", "pos"))
        
        #return sums
        return(data.table(perm=perm,
                          top=top,
                          sum=sum(gwas$Score.cline.beta),
                          sum..Z=sum(gwas$Score.cline.beta.Z.q),
                          sum.Z.=sum(gwas$Score.Z.q.cline.beta),
                          sum.Z.Z=sum(gwas$Score.Z.q.cline.beta.Z.q),
                          n=nrow(gwas),
                          n.conc=sum(gwas$conc, na.rm=T)))
        
    }
    return(rbindlist(a))
}

clump.enrich<-rbindlist(clump.enrich)
clump.enrich[,prop.conc:=n.conc/n]


write.table(clump.enrich, "/scratch/pae3g/evolution/index_snps_PRS_cline2019.txt", quote=F, sep="\t", row.names=F)

#now do 2014 cline


clump.enrich<-foreach(perm =c(0,101:200)) %do% { #cycle through original data and 100 permutations
    print(perm)
    a<-foreach (top = c(101, 490, 5316, 25250, 67914)) %do% { #these numbers correspond to the top # of snps in our different FDR cutoffs
        
        gwas<-fread(paste0("/scratch/pae3g/evolution/gwas_p_score_inv_id_perm", perm, ".txt"))
        snps<-fread(paste0("/nv/vol186/bergland-lab/Priscilla/perm_summaries/perm", perm, "_top", top, "_ld0.4_200kb.clumped"))
        cline<-fread("/scratch/pae3g/evolution/bergland_2014_cline_clean.txt")
        
        setnames(snps, c("CHR", "BP"), c("chr", "pos"))
        snps[,chr:=as.character(chr)]
        snps[chr=="23", chr:="X"]
        gwas<-merge(gwas, cline, by=c("chr", "pos"))
        #correct sign of cline for reference allele
        gwas[,clinal.beta:=clinal.beta*-1]
        gwas[,clinal.beta.q:=frank(clinal.beta)/(length(clinal.beta)+1)]
        gwas[,clinal.beta.Z.q:=qnorm(clinal.beta.q, 0,1)]
        
        #quantile rank and normalize Scores in merged datset
        gwas[,Score.q:=frank(gwas.Score)/(length(gwas.Score)+1)]
        gwas[,Score.Z.q:=qnorm(Score.q, 0,1)]
        #take product of Score * clinal beta
        gwas[,Score.cline.beta:=gwas.Score*clinal.beta]
        gwas[,Score.Z.q.cline.beta:=Score.Z.q*clinal.beta]
        gwas[,Score.Z.q.cline.beta.Z.q:=Score.Z.q*clinal.beta.Z.q]
        gwas[,Score.cline.beta.Z.q:=gwas.Score*clinal.beta.Z.q]
        
        
        #test for concordance of cline:
        gwas[,conc:=sign(clinal.beta)==sign(gwas.Score)]
        #take product of Score.Z*clinal.beta
        gwas<-merge(gwas, snps[, .(chr, pos)], by=c("chr", "pos"))
        
        #return sums
        return(data.table(perm=perm,
                          top=top,
                          sum=sum(gwas$Score.cline.beta),
                          sum..Z=sum(gwas$Score.cline.beta.Z.q),
                          sum.Z.=sum(gwas$Score.Z.q.cline.beta),
                          sum.Z.Z=sum(gwas$Score.Z.q.cline.beta.Z.q),
                          n=nrow(gwas),
                          n.conc=sum(gwas$conc, na.rm=T)))
        
    }
    return(rbindlist(a))
}

clump.enrich<-rbindlist(clump.enrich)
clump.enrich[,prop.conc:=n.conc/n]


write.table(clump.enrich, "/scratch/pae3g/evolution/index_snps_PRS_cline2014.txt", quote=F, sep="\t", row.names=F)


#now do seasonal stuff for machado 2018



clump.enrich<-foreach(perm =c(0,101:200)) %do% { #cycle through original data and 100 permutations
    print(perm)
    a<-foreach (top = c(101, 490, 5316, 25250, 67914)) %do% { #these numbers correspond to the top # of snps in our different FDR cutoffs
        
        gwas<-fread(paste0("/scratch/pae3g/evolution/gwas_p_score_inv_id_perm", perm, ".txt"))
        snps<-fread(paste0("/nv/vol186/bergland-lab/Priscilla/perm_summaries/perm", perm, "_top", top, "_ld0.4_200kb.clumped"))
        seas<-fread("/scratch/pae3g/evolution/seas_glm_switch_clean.txt")
        
        setnames(snps, c("CHR", "BP"), c("chr", "pos"))
        snps[,chr:=as.character(chr)]
        snps[chr=="23", chr:="X"]
        gwas<-merge(gwas, seas, by=c("chr", "pos"))
        #correct sign of cline for reference allele
        gwas[,seas.beta:=seas.beta*-1]
        gwas[,seas.beta.q:=frank(seas.beta)/(length(seas.beta)+1)]
        gwas[,seas.beta.Z.q:=qnorm(seas.beta.q, 0,1)]
        
        #quantile rank and normalize Scores in merged datset
        gwas[,Score.q:=frank(gwas.Score)/(length(gwas.Score)+1)]
        gwas[,Score.Z.q:=qnorm(Score.q, 0,1)]
        #take product of Score * clinal beta
        gwas[,Score.seas.beta:=gwas.Score*seas.beta]
        gwas[,Score.Z.q.seas.beta:=Score.Z.q*seas.beta]
        gwas[,Score.Z.q.seas.beta.Z.q:=Score.Z.q*seas.beta.Z.q]
        gwas[,Score.seas.beta.Z.q:=gwas.Score*seas.beta.Z.q]
        
        
        #test for concordance of cline:
        gwas[,conc:=sign(seas.beta)==sign(gwas.Score)]
        #take product of Score.Z*seas.beta
        gwas<-merge(gwas, snps[, .(chr, pos)], by=c("chr", "pos"))
        
        #return sums
        return(data.table(perm=perm,
                          top=top,
                          sum=sum(gwas$Score.seas.beta),
                          sum..Z=sum(gwas$Score.seas.beta.Z.q),
                          sum.Z.=sum(gwas$Score.Z.q.seas.beta),
                          sum.Z.Z=sum(gwas$Score.Z.q.seas.beta.Z.q),
                          n=nrow(gwas),
                          n.conc=sum(gwas$conc, na.rm=T)))
        
    }
    return(rbindlist(a))
}

clump.enrich<-rbindlist(clump.enrich)
clump.enrich[,prop.conc:=n.conc/n]


write.table(clump.enrich, "/scratch/pae3g/evolution/index_snps_PRS_seas2019.txt", quote=F, sep="\t", row.names=F)


##2014 seasonal data

clump.enrich<-foreach(perm =c(0,101:200)) %do% { #cycle through original data and 100 permutations
    print(perm)
    a<-foreach (top = c(101, 490, 5316, 25250, 67914)) %do% { #these numbers correspond to the top # of snps in our different FDR cutoffs
        
        gwas<-fread(paste0("/scratch/pae3g/evolution/gwas_p_score_inv_id_perm", perm, ".txt"))
        snps<-fread(paste0("/nv/vol186/bergland-lab/Priscilla/perm_summaries/perm", perm, "_top", top, "_ld0.4_200kb.clumped"))
        seas<-fread("/scratch/pae3g/evolution/bergland_2014_seasonal_clean.txt")
        
        setnames(snps, c("CHR", "BP"), c("chr", "pos"))
        snps[,chr:=as.character(chr)]
        snps[chr=="23", chr:="X"]
        gwas<-merge(gwas, seas, by=c("chr", "pos"))
        #no correciton in this dataset
        gwas[,seas.beta.q:=frank(seas.beta)/(length(seas.beta)+1)]
        gwas[,seas.beta.Z.q:=qnorm(seas.beta.q, 0,1)]
        
        #quantile rank and normalize Scores in merged datset
        gwas[,Score.q:=frank(gwas.Score)/(length(gwas.Score)+1)]
        gwas[,Score.Z.q:=qnorm(Score.q, 0,1)]
        #take product of Score * clinal beta
        gwas[,Score.seas.beta:=gwas.Score*seas.beta]
        gwas[,Score.Z.q.seas.beta:=Score.Z.q*seas.beta]
        gwas[,Score.Z.q.seas.beta.Z.q:=Score.Z.q*seas.beta.Z.q]
        gwas[,Score.seas.beta.Z.q:=gwas.Score*seas.beta.Z.q]
        
        
        #test for concordance of cline:
        gwas[,conc:=sign(seas.beta)==sign(gwas.Score)]
        #take product of Score.Z*seas.beta
        gwas<-merge(gwas, snps[, .(chr, pos)], by=c("chr", "pos"))
        
        #return sums
        return(data.table(perm=perm,
                          top=top,
                          sum=sum(gwas$Score.seas.beta),
                          sum..Z=sum(gwas$Score.seas.beta.Z.q),
                          sum.Z.=sum(gwas$Score.Z.q.seas.beta),
                          sum.Z.Z=sum(gwas$Score.Z.q.seas.beta.Z.q),
                          n=nrow(gwas),
                          n.conc=sum(gwas$conc, na.rm=T)))
        
    }
    return(rbindlist(a))
}

clump.enrich<-rbindlist(clump.enrich)
clump.enrich[,prop.conc:=n.conc/n]


write.table(clump.enrich, "/scratch/pae3g/evolution/index_snps_PRS_seas2014.txt", quote=F, sep="\t", row.names=F)


