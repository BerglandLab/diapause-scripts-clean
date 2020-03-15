#ihs analysis in dgrp around focal snps

#redownload dgrp

#start wtih alyssa's maine and pennsylvania files
cp /mnt/spicy_1/alyssa/inbredLines/vcf/BME_LN_LNPA.raw.vsqr.filter.recode.noREP.noINDEL.vcf /mnt/pricey_2/priscilla

#turn into phased haplotypes
sed -i 's/\//\|/g' /mnt/pricey_2/priscilla/BME_LN_LNPA.raw.vsqr.filter.recode.noREP.noINDEL.vcf

#get just genotypes
/mnt/pricey_2/priscilla/vcflib/bin/vcfkeepgeno /mnt/pricey_2/priscilla/BME_LN_LNPA.raw.vsqr.filter.recode.noREP.noINDEL.vcf GT > /mnt/pricey_2/priscilla/BME_LN_LNPA.gt.vcf
#split into body and header
grep -v "#" /mnt/pricey_2/priscilla/BME_LN_LNPA.gt.vcf > BME_LN_LNPA.vcf
grep "#" /mnt/pricey_2/priscilla/BME_LN_LNPA.gt.vcf > northern.header.txt

#replace missing data based on hardy weinberg (same script used for hybrid swarm)
swarm="north"
awk -v swarm=$swarm -v draw=$draw -v OFS="\t" 'BEGIN {
  print "chr\tpos\tnRef\tnAlt\tmostLikelyGeno\tn\ttotal">"/mnt/pricey_2/replace_bme_better.log"
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
{  printf $1"\t"$2"\t"nRef"\t"nAlt"\t"mostLikelyGeno"\t"n"\t"total"\n" > "/mnt/pricey_2/replace_bme_better.log"}
}' /mnt/pricey_2/priscilla/BME_LN_LNPA.vcf > /mnt/pricey_2/priscilla/BME_LN_LNPA.replaced.vcf

#concatenate back together
cat northern.header.txt /mnt/pricey_2/priscilla/BME_LN_LNPA.replaced.vcf > /mnt/pricey_2/priscilla/BME_LN_LNPA.replaced.final.vcf

#remove indels (there are alternate alleles with >1 nucleotide right now)

vcftools --vcf /mnt/pricey_2/priscilla/BME_LN_LNPA.replaced.final.vcf --remove-indels --recode --recode-INFO-all --out /mnt/pricey_2/priscilla/BME_LN_LNPA.replaced.final.noindel.vcf

#rename
mv BME_LN_LNPA.replaced.final.noindel.vcf.recode.vcf BME_LN_LNPA.replaced.final.noindel.vcf

grep -v "#" /mnt/pricey_2/priscilla/BME_LN_LNPA.replaced.final.noindel.vcf > /mnt/pricey_2/priscilla/BME_LN_LNPA.replaced.final.noindel.body.vcf
awk 'BEGIN { OFS = "\t"; } {$3 = $1"_"$2"_SNP"; print}'  /mnt/pricey_2/priscilla/BME_LN_LNPA.replaced.final.noindel.body.vcf >  /mnt/pricey_2/priscilla/BME_LN_LNPA.replaced.final.noindel.body.snpid.vcf

cat northern.header.txt /mnt/pricey_2/priscilla/BME_LN_LNPA.replaced.final.noindel.body.snpid.vcf > BME_LN_LNPA.final2.vcf

#use vcftools to make haplotype files
for chr in 2L 2R 3L 3R X; do
vcftools --vcf /mnt/pricey_2/priscilla/BME_LN_LNPA.final2.vcf --chr $chr --IMPUTE --max-missing 0.1 --out BME_LN_LNPA.replaced.final."$chr"
done

#remove headers from legends
#for chr in 2L 2R 3L 3R X; do
#  tail -n +2 dgrp2.filtered.impute."$chr".legend > dgrp2.filtered."$chr".inp
#done
