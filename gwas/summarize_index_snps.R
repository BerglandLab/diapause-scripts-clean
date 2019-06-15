library(data.table)
library(cowplot)
library(ggbeeswarm)
library(foreach)
clumps<-foreach(perm=c(0, 101:200))%do%{
    print(perm)
    #read in permutation p values and scores
    #loop through haplotype blocks
    e<-foreach(top=c(101, 490, 5316, 25250, 67914))%do%{
        blocks<-fread(paste0("/mnt/sammas_storage/bergland-lab/Priscilla/perm_summaries/perm", perm, "_top", top, "_ld0.4_200kb.clumped"))
        blocks[,perm:=perm]
        blocks[,top:=top]
       return(blocks)
    }
    e<-rbindlist(e)
}
clumps<-rbindlist(clumps)

write.table(clumps, "/mnt/pricey_2/priscilla/index_snp_perms.txt", sep="\t", quote=F, row.names=F)