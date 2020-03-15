bgzip -c hs.hc.A.keep.vcf > hs.hc.A.keep.vcf.gz
tabix -p vcf hs.hc.A.keep.vcf.gz

bgzip -c hs.hc.B.keep.vcf > hs.hc.B.keep.vcf.gz
tabix -p vcf hs.hc.B.keep.vcf.gz

bgzip -c hs.hc.snp.99.9.A.vcf > hs.hc.snp.99.9.A.vcf.gz
tabix -p vcf hs.hc.snp.99.9.A.vcf.gz

bgzip -c hs.hc.snp.99.9.B.vcf > hs.hc.snp.99.9.B.vcf.gz
tabix -p vcf hs.hc.snp.99.9.B.vcf.gz

vcf-merge hs.hc.A.keep.vcf.gz hs.hc.B.keep.vcf.gz > /mnt/pricey_2/priscilla/hybrid/hc/hs.hc.both.keep.vcf
vcf-merge hs.hc.snp.99.9.A.vcf.gz hs.hc.snp.99.9.B.vcf.gz > /mnt/pricey_2/priscilla/hybrid/hc/hs.hc.snp.99.9.both.vcf


#rearrange parent_line_by_swarm in R to reheader new vcfs

a<-fread("/mnt/pricey_2/priscilla/parent_line_by_swarm.txt", header=F)

a<-a[order(V3)][order(V4)]
write.table(a[,V1], "/mnt/pricey_2/priscilla/hybrid_swarm_line_names_A_B.txt", quote=F, row.names=F, sep="\t", col.names=F)

bcftools reheader -s /mnt/pricey_2/priscilla/hybrid_swarm_line_names_A_B.txt /mnt/pricey_2/priscilla/hybrid/hc/hs.hc.both.keep.vcf > //mnt/pricey_2/priscilla/hybrid/hc/hs.hc.both.keep.linenames.vcf

bcftools reheader -s /mnt/pricey_2/priscilla/hybrid_swarm_line_names_A_B.txt /mnt/pricey_2/priscilla/hybrid/hc/hs.hc.snp.99.9.both.vcf > //mnt/pricey_2/priscilla/hybrid/hc/hs.hc.snp.99.9.both.linenames.vcf

#find private alleles in parent vcfs

vcftools --vcf //mnt/pricey_2/priscilla/hybrid/hc/hs.hc.snp.99.9.both.linenames.vcf --out /mnt/pricey_2/priscilla/hybrid_swarm_parent --singletons
