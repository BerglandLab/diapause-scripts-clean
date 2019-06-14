#!/bin/bash
#
# #make covar file in R
# library(data.table)
# library(SNPRelate)
# phenos <- fread("/mnt/pricey_2/priscilla/hybrid/RivannaImputation/Genesis_files/phenotypes_with_calibration_031418.csv")
#  setkey(phenos, id)
#
#  ### some twiddles
#  phenos[,sex := "F"]
#  phenos[,id:=paste(Plate, Well, sep="_")]
#  phenos[,fid:=paste(Plate, Well, sep="_")]
#  phenos[,FID:=paste(Plate, Well, sep="_")]
#  phenos[,IID:=paste(Plate, Well, sep="_")]
#
#
# gcta=phenos[,.(id, fid, temp.rack.cal, photoperiod)]
# write.table(gcta, "/mnt/pricey_2/priscilla/qcov_gcta.txt", row.names=FALSE, quote=FALSE, col.names=FALSE, sep="\t")
# gcta2=phenos[,.(id, fid, swarm, generation)]
# write.table(gcta2, "/mnt/pricey_2/priscilla/cov_gcta.txt", row.names=FALSE, quote=FALSE, col.names=FALSE, sep="\t")
#  plink=phenos[,.(FID, IID, diapause.bin9)]
# write.table(plink, "/mnt/pricey_2/priscilla/plink_pheno.txt", sep="\t", quote=FALSE, row.names=FALSE)
#
# snpgdsVCF2GDS("/mnt/pricey_2/priscilla/final2.numchr2.vcf", "/mnt/pricey_2/priscilla/final2.numchr2.vcf.gds", method=c("biallelic.only"), snpfirstdim = FALSE)
# genofile <- snpgdsOpen("/mnt/pricey_2/priscilla/final2.numchr2.vcf.gds")
#
#     genotyped_samples<- read.gdsn(index.gdsn(genofile, "sample.id"))
#   #only use genotypes from swarm
#
#     ids.to.use=intersect(genotyped_samples, phenos$id)
#
# snps=data.table(snp.rs.id=read.gdsn(index.gdsn(genofile, "snp.rs.id")), snp.id=read.gdsn(index.gdsn(genofile, "snp.id")))
#
# snpset <- snpgdsLDpruning(genofile,
#                           ld.threshold=0.2,
#                           slide.max.bp = 5000,
#                           autosome.only=FALSE,
#                           missing.rate=.15,
#                           sample.id=ids.to.use,
#                           maf=.05)
# snpset<-unlist(snpset)
# snp<-snps[snp.id%in%snpset, snp.rs.id]
# write.table(snp, "/mnt/pricey_2/priscilla/ldsnpsforgatc.txt", quote=F, sep="\t", row.names=F, col.names=F)


### generate permuted phenotype, qcov, and covfiles for GATC


# library(data.table)
# library(foreach)
# phenos <- fread("/mnt/pricey_2/priscilla/hybrid/RivannaImputation/Genesis_files/phenotypes_with_calibration_031418.csv")
# setkey(phenos, id)
#
# perm.input<-fread("/mnt/pricey_2/priscilla/seeded_permutation_input.txt")
# foreach(perm=c(101:200))%do%{
# print(perm)
# seed<-perm.input[V4==1&V5==perm,V6]
# print(seed)
# set.seed(seed)
# phenos[swarm=="A", id:=sample(id)]
# phenos[swarm=="B", id:=sample(id)]
# write.table(phenos[,.(id, id, diapause.bin9)], paste0("/mnt/pricey_2/priscilla/plink_phenos_perm", perm, ".txt"), quote=F, sep="\t", row.names=F, col.names=F)
# write.table(phenos[,.(id, id, swarm, generation)], paste0("/mnt/pricey_2/priscilla/cov_gcta__perm", perm, ".txt"), quote=F, sep="\t", row.names=F, col.names=F)
# write.table(phenos[,.(id, id, temp.rack.cal, photoperiod)], paste0("/mnt/pricey_2/priscilla/qcov_gcta__perm", perm, ".txt"), quote=F, sep="\t", row.names=F, col.names=F)
# write.table(phenos, paste0("/mnt/pricey_2/priscilla/phenos_perm", perm, ".txt"), quote=F, sep="\t", row.names=F, col.names=T)
#  }
# #
# ##need to renumber chromosomes in vcf
#



cp /mnt/pricey_2/priscilla/final2.vcf /mnt/pricey_2/priscilla/final2.numchr.vcf

sed -i 's/X/1/g' /mnt/pricey_2/priscilla/final2.numchr.vcf
sed -i 's/2L/2/g' /mnt/pricey_2/priscilla/final2.numchr.vcf
sed -i 's/2R/3/g' /mnt/pricey_2/priscilla/final2.numchr.vcf
sed -i 's/3L/4/g' /mnt/pricey_2/priscilla/final2.numchr.vcf
sed -i 's/3R/5/g' /mnt/pricey_2/priscilla/final2.numchr.vcf

grep -v -e Yhet -e YHet -e gi -e dmel /mnt/pricey_2/priscilla/final2.numchr.vcf > /mnt/pricey_2/priscilla/final2.numchr2.vcf


plink1.9 --vcf /mnt/pricey_2/priscilla/final2.numchr2.vcf --double-id --vcf-half-call missing --out /mnt/pricey_2/priscilla/final2.numchr2


#make GRM
/mnt/pricey_2/priscilla/gcta_1.91.3beta/gcta64 --bfile /mnt/pricey_2/priscilla/final2.numchr2 --make-grm --out /mnt/pricey_2/priscilla/final2.grm --thread-num 10

#make GRM with LDpruned snps
/mnt/pricey_2/priscilla/gcta_1.91.3beta/gcta64 --bfile /mnt/pricey_2/priscilla/final2.numchr2 --make-grm --out /mnt/pricey_2/priscilla/final2.snpset.grm --thread-num 10 --extract /mnt/pricey_2/priscilla/ldsnpsforgatc.txt


#calculate heritability
for perm in {101..200}; do
  echo $perm
  /mnt/pricey_2/priscilla/gcta_1.91.3beta/gcta64 \
  --grm /mnt/pricey_2/priscilla/final2.snpset.grm \
  --pheno /mnt/pricey_2/priscilla/plink_phenos_perm$perm.txt \
  --mpheno 1 --qcovar /mnt/pricey_2/priscilla/qcov_gcta__perm$perm.txt \
  --covar /mnt/pricey_2/priscilla/cov_gcta__perm$perm.txt \
  --reml \
  --reml-alg 1 \
  --reml-bendV \
  --out /mnt/pricey_2/priscilla/gcta_$perm \
  --thread-num 10
done

#exctract heritability in R
library(data.table)
library(foreach)

y<-foreach(perm=c(101:200))%do%{
  h<-fread(paste0("/mnt/pricey_2/priscilla/gcta_", perm, ".hsq"), fill=T)
  h[,perm:=perm]
  return(h)
}

y<-rbindlist(y)

a<-fread("/mnt/pricey_2/priscilla/gcta.hsq", fill=T)
a[,perm:=0]

y<-rbind(a, y)

write.table(y,"/mnt/pricey_2/priscilla/gcta_hsq.txt", quote=F, sep="\t", row.names=F)

ggplot(y[Source=="V(G)/Vp"&perm!=0], aes(x=Variance))+geom_histogram(fill="grey50", binwidth=0.001)+geom_vline(xintercept=a[Source=="V(G)/Vp", Variance], color="lightseagreen", size=2)+labs(x="VG/VP")
