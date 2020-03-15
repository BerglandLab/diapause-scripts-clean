#annotate GWAS SNPs using snpEff

#on workstation in /mnt/pricey_2/priscilla
wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
unzip snpEff_latest_core.zip

# want to use the database called BDGP5.75

java -Xmx4g -jar /mnt/pricey_2/priscilla/snpEff/snpEff.jar BDGP5.75 /mnt/pricey_2/priscilla/all_sites.vcf > all_sites.ann.vcf

java -Xmx4g -jar /mnt/pricey_2/priscilla/snpEff/snpEff.jar \
closest BDGP5.75 /mnt/pricey_2/priscilla/all_sites.vcf \
  -canon \
  -t \
  > all_sites.ann.closest.canon.vcf


  java -Xmx4g -jar /mnt/pricey_2/priscilla/snpEff/snpEff.jar \
   BDGP5.75 /mnt/pricey_2/priscilla/all_sites.vcf \
    -canon \
    -t \
    > all_sites.ann.canon.vcf
