#install.packages("rehh")
# library(rehh)
# library(data.table)
# library(foreach)
# library(doMC)
# registerDoMC(20)
#fix input files
# foreach(chr=c("2L", "2R", "3L", "3R", "X"))%do%{
#     inp<-fread(paste0("/mnt/pricey_2/priscilla/dgrp2.filtered.", chr, ".impute.legend"), header=T)
#     inp[,chr:=tstrsplit(ID, split="_")[[1]]]
#     inp[,allele0:=0]
#     inp[,allele1:=1]
#     write.table(inp[,.(ID, chr, pos, allele0, allele1)], paste0("/mnt/pricey_2/priscilla/dgrp2.", chr, ".inp"), quote=F, sep=" ", row.names=F, col.names=F)
# }
# 
# a<-foreach(chr=c("2L", "2R", "3L", "3R", "X"))%do%{
#     hap<-data2haplohh(hap_file=paste0("dgrp2.filtered.", chr,".impute.hap"),map_file=paste0("dgrp2.", chr, ".inp"),min_perc_geno.hap=90, min_perc_geno.snp=90, recode.allele=TRUE, haplotype.in.columns=TRUE)
# 
#     res<-scan_hh(hap, threads=20)
#     return(as.data.table(res))
# }
# a<-rbindlist(a)
# 
# write.table(a, "/mnt/pricey_2/priscilla/ihs.txt", quote=F, sep="\t", row.names=F)

#move ihs.txt to rivanna to work there
#flip polarity according to pro-diapause allele

library(data.table)
library(foreach)
library(rehh)
 library(doMC)
registerDoMC(20)

a.raw<-fread("/mnt/pricey_2/priscilla/ihs.txt")

th.test<-foreach(perm=c(0, 101:200))%do%{
    print(perm)
    #read in permutation p values and scores
    gwas<-fread(paste0("/mnt/sammas_storage/bergland-lab/Priscilla/perm_summaries/gwas_p_score_inv_id_perm", perm, ".txt"))


    setnames(gwas, c("chr", "pos"), c("CHR", "POSITION"))

    a<-merge(a.raw, gwas, by=c("CHR", "POSITION"), all.x=T)
#currently coded so that ancestral allele is reference (0) and derived allele is alternate (1)
#we want "derived" to be the pro-diapause allele
#positive score means that ref is pro-diapause (aka ancestral is pro diapause)
#so, when gwas score is POSTIVE, switch ihh_a and ihh_d before computing ihs

    a[gwas.Score>=0,ihh.a:=iHH_D]
    a[is.na(gwas.Score)|gwas.Score<0,ihh.a:=iHH_A]
    a[gwas.Score>=0,ihh.d:=iHH_A]
    a[is.na(gwas.Score)|gwas.Score<0,ihh.d:=iHH_D]

#flip allele frequency if ancestral and derived are flipped
    a[gwas.Score>=0, freq_A:=1-freq_A]

#delete original calls and rename
    a[,iHH_A:=NULL]
    a[,iHH_D:=NULL]
    setnames(a, c("ihh.a", "ihh.d"), c("iHH_A", "iHH_D"))
    a<-a[,.(CHR, POSITION, freq_A,  iHH_A, iHH_D, iES_Tang_et_al_2007, iES_Sabeti_et_al_2007)]

    b<-ihh2ihs(a, minmaf=0)
    c<-as.data.table(b$iHS)
    
    d<-merge(gwas, c, by=c("CHR", "POSITION"))
    setnames(d, c("CHR","POSITION", "-log10(p-value)"), c("chr", "pos",'ihs.p'))


    e<-foreach(top=c(101, 490, 5316, 25250, 67914))%dopar%{
        blocks<-fread(paste0("/mnt/sammas_storage/bergland-lab/Priscilla/perm_summaries/perm", perm, "_top", top, "_ld0.4_200kb.clumped"))
        setnames(blocks, c("CHR", "BP"), c("chr", "pos"))
        blocks[,chr:=as.character(chr)]
        blocks[chr=="23", chr:="X"]
        f<-merge(blocks[,.(chr, pos)], d, by=c("chr", "pos"))
        #pull out most significant snp in each block
        f[,perm:=perm]
        f[,top:=top]
        return(f)
    }

    e<-rbindlist(e)
    return(e)
}

th.test<-rbindlist(th.test)
write.table(th.test, "/mnt/pricey_2/priscilla/ihs_dgrp_clump200kb.txt", quote=F, row.names=F, sep="\t")

library(data.table)
library(ggbeeswarm)
library(cowplot)
th.test<-fread("/mnt/pricey_2/priscilla/ihs_dgrp_clump200kb.txt")

th.test.sum<-th.test[,
                     .(n=.N,
                       min.ihs=min(iHS, na.rm=T), 
                        med.ihs=median(iHS, na.rm=T), 
                        max.ihs=max(iHS, na.rm=T)), 
                     .(top, perm)]

th.test.sum.melt<-melt(th.test.sum, id.vars=c("top", "perm", "n"))

ggplot(th.test.sum.melt[perm!=0], aes(x=as.factor(top), y=value, color=variable))+
    geom_quasirandom(dodge.width = .8, method="smiley")+
    geom_point(data=th.test.sum.melt[perm==0], aes(x=as.factor(top), y=value, group=variable), color="black" ,position=position_dodge(width=0.8))+
    labs(x="Top # SNPs", y="iHS", color="iHS of most significant\nGWAS SNP in haplotype")+
    scale_color_discrete(labels=c("minimum", "median", "maximum"))

