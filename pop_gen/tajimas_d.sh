#tajima's d
grep "#" dgrp2.vcf > header.txt
cat header.txt dgrp2.filtered.vcf > dgrp2.fh.vcf
#make files and put in folder in bash:

 for chr in 2L 2R 3L 3R X; do
  echo $chr
  cd /mnt/pricey_2/priscilla/
  rm -r $chr
  grep -e "#" -e $chr dgrp2.fh.vcf > dgrp2.$chr.vcf
  mkdir $chr
  mv dgrp2.$chr.vcf $chr
  cd $chr
  bgzip dgrp2.$chr.vcf
  tabix -p vcf dgrp2.$chr.vcf.gz
 done
