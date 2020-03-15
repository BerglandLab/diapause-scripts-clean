#!/bin/bash


#make GRM using all snps
/mnt/pricey_2/priscilla/gcta/gcta_1.91.3beta/gcta64 --bfile /mnt/pricey_2/priscilla/final2.numchr2 --make-grm --out /mnt/pricey_2/priscilla/final2.grm --thread-num 10

#calculate heritability
for perm in {1000..1999}; do
  echo $perm
  /mnt/pricey_2/priscilla/gcta/gcta_1.91.3beta/gcta64 \
  --grm /mnt/pricey_2/priscilla/final2.grm \
  --pheno /mnt/pricey_2/priscilla/gcta/plink_phenos_perm$perm.txt \
  --mpheno 1 --qcovar /mnt/pricey_2/priscilla/gcta/qcov_gcta_perm$perm.txt \
  --covar /mnt/pricey_2/priscilla/gcta/cov_gcta_perm$perm.txt \
  --reml \
  --reml-alg 1 \
  --reml-bendV \
  --out /mnt/pricey_2/priscilla/gcta/gcta_$perm \
  --thread-num 10
done

for perm in {1000..1999}; do
  echo $perm
  /mnt/pricey_2/priscilla/gcta/gcta_1.91.3beta/gcta64 \
  --grm /mnt/pricey_2/priscilla/final2.grm \
  --pheno /mnt/pricey_2/priscilla/gcta/plink_phenos_perm$perm.txt \
  --mpheno 2 --qcovar /mnt/pricey_2/priscilla/gcta/qcov_gcta_perm$perm.txt \
  --covar /mnt/pricey_2/priscilla/gcta/cov_gcta_perm$perm.txt \
  --reml \
  --reml-alg 1 \
  --reml-bendV \
  --out /mnt/pricey_2/priscilla/gcta/gcta_st7$perm \
  --thread-num 10
done
