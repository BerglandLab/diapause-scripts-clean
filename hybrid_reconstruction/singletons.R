
library(data.table)
s<-fread("/mnt/pricey_2/priscilla/hybrid_swarm_parent.singletons")
f<-fread("/mnt/pricey_2/priscilla/hwe_missing_maf_filters.txt")
setnames(s, c("CHROM", "POS"), c("chr", "pos"))

x<-merge(f, s, by=c("chr", "pos"))

x[freq==1]

#only 109 of 806,000 singletons were lost in the hybrid swarms. need to work this in somewhere.
s[chr%in%c("2L", "3L", "2R", "3R", "X")]
