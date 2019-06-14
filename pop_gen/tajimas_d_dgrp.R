#tajima's d
#grep "#" dgrp2.vcf > header.txt
#cat header.txt dgrp2.filtered.vcf > dgrp2.fh.vcf
#make files and put in folder in bash:
# 
# for chr in 2L 2R 3L 3R X; do
# echo $chr
# cd /mnt/pricey_2/priscilla/
# rm -r $chr
# grep -e "#" -e $chr dgrp2.fh.vcf > dgrp2.$chr.vcf
# mkdir $chr
# mv dgrp2.$chr.vcf $chr
# cd $chr
# bgzip dgrp2.$chr.vcf
# tabix -p vcf dgrp2.$chr.vcf.gz
# done

#install.packages("PopGenome")
# library(PopGenome)
# library(foreach)
# library(data.table)
# 
# 
# dgrpstats<-foreach(chr=c("2L", "2R", "3L", "3R", "X"))%do%{
#     vcf<-readVCF(paste0("/mnt/pricey_2/priscilla/", chr, "/dgrp2.", chr, ".vcf.gz"), numcols=10000, tid=chr, frompos=1, topos=30000000, include.unknown=TRUE)
#     sliding <- sliding.window.transform(vcf,1000,500,type=2)
#     vcf <- neutrality.stats(vcf, FAST=T)
#     sliding <- neutrality.stats(sliding)
#     neutrality <- get.neutrality(sliding)[[1]]
#     sliding <- linkage.stats(sliding, do.ZnS = T)
#     linkage <- get.linkage(sliding)[[1]]
#     sliding <- diversity.stats(sliding)
#     diversity <- get.diversity(sliding)[[1]]
#     stats <- cbind(neutrality,linkage,diversity)
#     #clean up stats dataframe
#     stats <- cbind(Row.names=rownames(stats),stats)
#     rownames(stats) <- NULL
#     stats <- as.data.table(stats)
#     stats[,Start:=as.numeric(tstrsplit(Row.names, split=" - ")[[1]])]
#     stats[,End:=Start+999]
#     stats[,CHROM:=chr]
#     rm(vcf)
#     return(stats)
# }
# 
# dgrpstats<-rbindlist(dgrpstats)
# 
# write.table(dgrpstats, "/mnt/pricey_2/priscilla/dgrp_tajimasd.txt", quote=F, row.names=F, sep="\t")


#analyze with haplotype blocks
library(PopGenome)
library(foreach)
library(data.table)

dgrpstats<-fread("/mnt/pricey_2/priscilla/dgrp_tajimasd.txt")

#make a table with just TJ statistic and chromosome info
tj<-dgrpstats[,.(CHROM, Start, End, Tajima.D)]
setnames(tj, "CHROM", "chr")

#loop through permutations
th.test<-foreach(perm=c(0, 101:200))%do%{
    print(perm)
    #read in permutation p values and scores
    gwas<-fread(paste0("/mnt/sammas_storage/bergland-lab/Priscilla/perm_summaries/gwas_p_score_inv_id_perm", perm, ".txt"))
    #loop through haplotype blocks
    e<-foreach(top=c(101, 490, 5316, 25250, 67914))%do%{
        blocks<-fread(paste0("/mnt/sammas_storage/bergland-lab/Priscilla/perm_summaries/perm", perm, "_top", top, "_ld0.4_200kb.clumped"))
        setnames(blocks, c("CHR", "BP"), c("chr", "pos"))
        blocks[,chr:=as.character(chr)]
        blocks[chr=="23", chr:="X"]
        f<-merge(blocks, gwas, by=c("chr", "pos"))
        #pull out most significant snp in each block
        #get the tajima's d values for the most significant snps (2 values per snp because every snp appears in two overlapping windows)
        tajima.test<-foreach(i=c(1:nrow(f)))%do%{
            test.chr=f[i,chr]
            test.pos=f[i,pos]
            tajimas=tj[chr==test.chr&Start<test.pos&End>test.pos, Tajima.D]
            return(data.table(i=i,
                              chr=test.chr,
                              pos=test.pos,
                              tajimas.d=tajimas))
        }
        tajima.test<-rbindlist(tajima.test)
        tajima.test[,perm:=perm]
        tajima.test[,top:=top]
        return(tajima.test)
    }
    e<-rbindlist(e)
}
th.test<-rbindlist(th.test)

write.table(th.test, "/mnt/pricey_2/priscilla/tajimasd_clumps200kb.txt", quote=F, row.names=F, sep="\t")

tj<fread("/mnt/pricey_2/priscilla/tajimasd_clumps200kb.txt")

tj.sum<-tj[,.(med=median(as.numeric(tajimas.d), na.rm=T),
                   max=max(as.numeric(tajimas.d), na.rm=T),
                   min=min(as.numeric(tajimas.d), na.rm=T)), 
                .(perm, top)]
tj.sum.melt<-melt(tj.sum, id.vars=c("perm", "top"))

ggplot(tj.sum.melt[perm!=0], aes(x=as.factor(top), y=value, color=variable))+
    geom_quasirandom(dodge.width = .8, method="smiley")+
    geom_point(data=tj.sum.melt[perm==0], aes(x=as.factor(top), y=value, group=variable), color="black" ,position=position_dodge(width=0.8))+
    labs(x="Top # SNPs", y="tajima's d", color="")#+
    #scale_color_discrete(labels=c("minimum", "median", "maximum"))
