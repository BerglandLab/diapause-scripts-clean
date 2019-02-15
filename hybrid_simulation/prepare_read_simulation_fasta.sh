#!/bin/bash

#first split most conservative vcf into individual files in parallel using vcf-subset

bcftools query -l /mnt/pricey_2/priscilla/hybrid/hc/hs.hc.A.keep.vcf | parallel -j 8 'vcf-subset -c {} /mnt/pricey_2/priscilla/hybrid/hc/hs.hc.A.keep.vcf > /mnt/pricey_2/priscilla/hybrid/hc/simulation/{}.A.keep.vcf'

bcftools query -l /mnt/pricey_2/priscilla/hybrid/hc/hs.hc.B.keep.vcf | parallel -j 8 'vcf-subset -c {} /mnt/pricey_2/priscilla/hybrid/hc/hs.hc.B.keep.vcf > /mnt/pricey_2/priscilla/hybrid/hc/simulation/{}.B.keep.vcf'

ls  /mnt/pricey_2/priscilla/hybrid/hc/simulation/*.vcf >  /mnt/pricey_2/priscilla/hybrid/hc/simulation/single_vcfs.txt

while read file ; do
    #echo $file
    filename=${file##*/}
    #echo $filename
    file_id=${filename%.vcf}
    echo $file_id
    sample=${filename%.*.*.vcf}
    #echo $sample
    for chr in 2L 2R 3L 3R X; do
	echo $chr
	java -jar /usr/local/bin/GenomeAnalysisTK.jar \
	     -T FastaAlternateReferenceMaker \
	     -R /mnt/pricey_2/priscilla/hybrid/seqs/all_dmel.fasta   \
	     -o /mnt/pricey_2/priscilla/hybrid/hc/simulation/$file_id.$chr.fasta \
	     -V $file \
	     -L $chr \
	     --use_IUPAC_sample $sample \
	     --rawOnelineSeq 
    done 
done < /mnt/pricey_2/priscilla/hybrid/hc/simulation/single_vcfs.txt


#rename these files and move into separate folders sample_Chr.seq

while read line; do 
    for swarm in A B; do
	for chr in 2L 2R 3L 3R X; do
	    mkdir -p /mnt/pricey_2/priscilla/hybrid/hc/simulation/$chr/$swarm
	    cp /mnt/pricey_2/priscilla/hybrid/hc/simulation/$line.$swarm.keep.$chr.fasta /mnt/pricey_2/priscilla/hybrid/hc/simulation/$chr/$swarm/$line"_Chr"$chr.seq
	done
    done
done < /mnt/pricey_2/priscilla/hybrid/etc/34lines.txt




