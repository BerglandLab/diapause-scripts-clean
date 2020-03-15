#ihs analysis in dgrp around focal snps

#redownload dgrp

wget http://dgrp2.gnets.ncsu.edu/data/website/dgrp2.
#turn into phased haplotypes
sed -i 's/\//\|/g' /mnt/pricey_2/priscilla/dgrp2.vcf
#split into body and header
grep SNP dgrp2.vcf > dgrp2.filtered.vcf
grep "#" dgrp2.vcf > header.txt

#replace missing data based on hardy weinberg (same script used for hybrid swarm)
swarm="dgrp2"
awk -v swarm=$swarm -v draw=$draw -v OFS="\t" 'BEGIN {
  print "chr\tpos\tnRef\tnAlt\tmostLikelyGeno\tn\ttotal">"/mnt/pricey_2/replace_dgrp2_better.log"
}
{
  nRef=0
  nAlt=0
  for(i=1; i<=NF; i++) {
if($i=="0|0") nRef+=2
if($i=="0|." || $i==".|0") nRef+=1
if($i=="1|." || $i==".|1") nAlt+=1
if($i=="0|1" || $i=="1|0"){
    nRef+=1
    nAlt+=1
}
if($i=="1|1") nAlt+=2
  }
  if((nRef+nAlt)==0){
refFreq=1
altFreq=0
  }
  else{
refFreq=nRef/(nRef+nAlt)
altFreq=1-refFreq
  }
  freq00=refFreq*refFreq
  freq01=refFreq*altFreq*2
  freq11=altFreq*altFreq
  if(freq00 > freq01 && freq00 > freq11){
mostLikelyGeno="homRef"
n=gsub(/\./, 0); print
}
  if(freq01 > freq00 && freq01 > freq11){
mostLikelyGeno="het"
n=gsub(/0\|\./, "0|1") + gsub(/\.\|0/, "1|0") + gsub(/1\|\./, "1|0") + gsub(/\.\|1/, "0|1")+gsub(/\.\|\./, "0|1")*2; print
  }
  if(freq11 > freq01 && freq11>freq00){
mostLikelyGeno="homAlt"
n=gsub(/\./, 1); print
  }
total=nRef+nAlt+n
{  printf $1"\t"$2"\t"nRef"\t"nAlt"\t"mostLikelyGeno"\t"n"\t"total"\n" > "/mnt/pricey_2/replace_dgrp2_better.log"}
}' /mnt/pricey_2/priscilla/dgrp2.filtered.vcf > /mnt/pricey_2/priscilla/dgrp2.filtered.replaced.vcf

#concatenate back together
cat header.txt dgrp2.filtered.replaced.vcf > dgrp2.filtered.final.txt

#use vcftools to make haplotype files
for chr in 2L 2R 3L 3R X; do
vcftools --vcf dgrp2.filtered.final.txt --chr $chr --IMPUTE --max-missing 0.1 --out dgrp2.filtered."$chr"
done

#remove headers from legends
#for chr in 2L 2R 3L 3R X; do
#  tail -n +2 dgrp2.filtered.impute."$chr".legend > dgrp2.filtered."$chr".inp
#done
