#!/bin/bash

#starting over with haplotype caller to call snps and indels, then will clean up 

#for haplotype caller need to have a list of bams for each gvcf. this is messy because the files have all sorts of names, but can grep through all possiblenames to get each

#have to go to rivanna because gvcf creation is super slow and memory-intensive on workstation.



####
#HAPLOTYPE CALLER ON RIVANNA #
#######

#parallelized rsync to move bams

cat /mnt/internal_2/priscilla/HSparents/scripts/mergedbams.list | parallel --will-cite -j 12 rsync -avzhe ssh {} pae3g@rivanna1.hpc.virginia.edu:/scratch/pae3g/parents/bams

#also need to move all bam.bai files over to rivanna. first add ".bai" to the end of every file
rsync -avzhe ssh /mnt/internal_2/priscilla/HSparents/mapped/*.bai pae3g@rivanna1.hpc.virginia.edu:/scratch/pae3g/parents/bams

#copy bam names and file names to rivanna
scp /mnt/internal_2/priscilla/HSparents/inbred_mapped/inbred_all_names.txt pae3g@rivanna1.hpc.virginia.edu:/scratch/pae3g/parents/bams
scp /mnt/internal_2/priscilla/HSparents/inbred_mapped/inbred_new_names.txt pae3g@rivanna1.hpc.virginia.edu:/scratch/pae3g/parents/bams
	
#on rivanna, make a list of all bams:
ls /scratch/pae3g/parents/bams/*.bam > /scratch/pae3g/parents/bams/mergedbams.list

#make file listing all bams for each file using both old and new file names
while read old new filename; do
    echo 'old name =' $old
    echo 'new name =' $new
    grep -e $old -e $new -e $filename /scratch/pae3g/parents/bams/mergedbams.list > /scratch/pae3g/parents/bams/"$new".bams.list 
done < /scratch/pae3g/parents/bams/inbred_all_names.txt

#loop through batch submit to do haplotype caller
while read line; do
      sbatch /scratch/pae3g/scripts/submit_haplotype_caller.slurm "$line"
done </scratch/pae3g/parents/bams/inbred_new_names.txt

#move gvcfs back to workstation (from workstation)

nohup rsync -avzh pae3g@rivanna1.hpc.virginia.edu:/scratch/pae3g/parents/gVCF /mnt/internal_2/priscilla/HSparents/variants/ &

#DGRP 28243 is going to take forever on workstation. make multi-threaded on workstation

####
###gVCF ASSEMBLY ON WORKSTATION####
####

#downsample this file
samtools view -s .2 -b DGRP_28243.SRR835347.sort.dedup.renamed.bam > DGRP_28243.SRR835347.sort.dedup.renamed.downsampled.bam

#run haplotype caller on downsampled file
nohup java -jar /usr/local/bin/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R /mnt/internal_2/priscilla/ref_files/all_dmel.fasta \
-I /mnt/internal_2/priscilla/HSparents/scripts/DGRP_28243edit.bams.list \
-hets 0.01 \
-indelHeterozygosity 0.001 \
--emitRefConfidence GVCF \
--variant_index_type LINEAR \
--variant_index_parameter 128000 \
-o /mnt/internal_2/priscilla/HSparents/variants/gVCF/DGRP_28243downsampled.raw.snps.indels.g.vcf &


#make list of all gvcf files
ls /mnt/internal_2/priscilla/HSparents/variants/gVCF/*.g.vcf > /mnt/internal_2/priscilla/HSparents/variants/gVCF/gvcfs.list

paste <(sort -k1,1 /mnt/internal_2/priscilla/HSparents/inbred_mapped/inbred_new_names.txt) <(sort -k1,1 /mnt/internal_2/priscilla/HSparents/variants/gVCF/gvcfs.list) > /mnt/internal_2/priscilla/HSparents/variants/gVCF/merged.hc_map

#combine into one vcf

#install gatk4 locally

#git clone https://github.com/broadinstitute/gatk.git

#export PATH=$PATH:/mnt/internal_2/priscilla/HSparents/gatk-4.0.1.2/gatk

#for chr in 2L 2R 3L 3R X; do

#/mnt/internal_2/priscilla/HSparents/gatk-4.0.1.2/gatk --java-options "-Xmx4G" \
#       GenomicsDBImport \
#       --genomicsdb-workspace-path /mnt/internal_2/priscilla/HSparents/genomicsDB_"$chr" \
#       --sample-name-map /mnt/internal_2/priscilla/HSparents/variants/gVCF/merged.hc_map \
#       --L $chr 
 #      done


#genotype gvcfs from first round of genomicsDB
#/mnt/internal_2/priscilla/HSparents/gatk-4.0.1.2/gatk --java-options "-Xmx4G" \
#     GenotypeGVCFs \
#     -R /mnt/internal_2/priscilla/ref_files/all_dmel.fasta \
#     -V gendb:///mnt/internal_2/priscilla/HSparents/genomicsDB_2L \
#     -O /mnt/internal_2/priscilla/HSparents/variants/merged.2L.gdb.hc.raw.vcf


#combine all vcfs with CombineGVCFs  and genotype(not using genomicsDB because it didn't complete)
nohup java -jar /usr/local/bin/GenomeAnalysisTK.jar \
      -T CombineGVCFs \
      -R /mnt/internal_2/priscilla/ref_files/all_dmel.fasta \
      -V /mnt/internal_2/priscilla/HSparents/variants/gVCF/gvcfs.list \
      -o /mnt/internal_2/priscilla/HSparents/variants/gVCF/merged.hc.combined.g.vcf &


#genotype (lines of raw file = 4949593)
nohup java -jar /usr/local/bin/GenomeAnalysisTK.jar \
     -T GenotypeGVCFs \
     -R /mnt/internal_2/priscilla/ref_files/all_dmel.fasta \
     -nt 4 \
     -V /mnt/internal_2/priscilla/HSparents/variants/gVCF/merged.hc.combined.g.vcf \
     -o /mnt/internal_2/priscilla/HSparents/variants/merged.hc.raw.vcf &

### make thruthiness set using Alan's random dgrp subset file
bedtools intersect -header \
	-a /mnt/internal_2/priscilla/HSparents/variants/merged.hc.raw.vcf \
	-b /mnt/internal_2/priscilla/HSparents/variants/dgrp2_snp.randomSubset.noREP.noINDEL.vcf > \
	/mnt/internal_2/priscilla/HSparents/variants/merged.hc.dgrp2_snp.randomSubset.noREP.noINDEL.vcf.truthiness.vcf

java -jar /usr/local/bin/GenomeAnalysisTK.jar \
      -T VariantRecalibrator \
      -R /mnt/internal_2/priscilla/ref_files/all_dmel.fasta \
      -input /mnt/internal_2/priscilla/HSparents/variants/merged.hc.raw.vcf \
      -resource:dgrp,known=false,training=true,truth=true,prior=15.0 /mnt/internal_2/priscilla/HSparents/variants/merged.hc.dgrp2_snp.randomSubset.noREP.noINDEL.vcf.truthiness.vcf \
      -an DP \
      -an QD \
      -an FS \
      -an SOR \
      -an MQ \
      -an MQRankSum \
      -an ReadPosRankSum \
      -mode SNP \
      -nt 6 \
      -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
      -recalFile /mnt/internal_2/priscilla/HSparents/variants/hc.recalibrate_SNP_merged.recal \
      -tranchesFile /mnt/internal_2/priscilla/HSparents/variants/hc.recalibrate_SNP.tranches \
      -rscriptFile /mnt/internal_2/priscilla/HSparents/variants/hc.recalibrate_SNP_plots_merged.R 

#apply recalibration (4949601)
java -jar /usr/local/bin/GenomeAnalysisTK.jar \
      -T ApplyRecalibration \
      -R /mnt/internal_2/priscilla/ref_files/all_dmel.fasta \
      -input /mnt/internal_2/priscilla/HSparents/variants/merged.hc.raw.vcf \
      -mode SNP \
      --ts_filter_level 99.0 \
      -recalFile /mnt/internal_2/priscilla/HSparents/variants/hc.recalibrate_SNP_merged.recal \
      -tranchesFile /mnt/internal_2/priscilla/HSparents/variants/hc.recalibrate_SNP.tranches \
      -o /mnt/internal_2/priscilla/HSparents/variants/merged.hc.raw.vsqr.vcf

#get file of all snps that did not pass (n=818687)

grep VQSRTrancheSNP /mnt/internal_2/priscilla/HSparents/variants/merged.hc.raw.vsqr.vcf > /mnt/internal_2/priscilla/HSparents/variants/merged.hc.raw.vsqr.tranch99.removedSNPs

for chr in 2L 2R 3L 3R X ; do
    echo $chr
    grep $chr /mnt/internal_2/priscilla/HSparents/variants/merged.hc.raw.vsqr.tranch99.removedSNPs | wc -l
done


#2L
#177429
#2R
#176909
#3L
#200758
#3R
#132136
#X
#102924

#clean up with filter (4130917 variants after, still includes indels and mnps)
vcftools \
--vcf /mnt/internal_2/priscilla/HSparents/variants/merged.hc.raw.vsqr.vcf \
--remove-filtered-all \
--recode \
--recode-INFO-all \
--out /mnt/internal_2/priscilla/HSparents/variants/merged.hc.filter.vsqr

########
# MORE LENIENT FILTRATION POSSIBILITY
#######

#apply recalibration at a more lenient threshold to increase snp count (4949600)
java -jar /usr/local/bin/GenomeAnalysisTK.jar \
      -T ApplyRecalibration \
      -R /mnt/internal_2/priscilla/ref_files/all_dmel.fasta \
      -input /mnt/internal_2/priscilla/HSparents/variants/merged.hc.raw.vcf \
      -mode SNP \
      --ts_filter_level 99.9 \
      -recalFile /mnt/internal_2/priscilla/HSparents/variants/hc.recalibrate_SNP_merged.recal \
      -tranchesFile /mnt/internal_2/priscilla/HSparents/variants/hc.recalibrate_SNP.tranches \
      -o /mnt/internal_2/priscilla/HSparents/variants/merged.hc.raw.vsqr.99.9.vcf  



#clean up with filter with more lenient 99.9% tranch( 4410943 variants after, so this keeps extra ~300K snps)
vcftools \
--vcf /mnt/internal_2/priscilla/HSparents/variants/merged.hc.raw.vsqr.99.9.vcf \
--remove-filtered-all \
--recode \
--recode-INFO-all \
--out /mnt/internal_2/priscilla/HSparents/variants/merged.hc.filter.99.9.vsqr

#count removed snps (538540 removed at more lenient threshold)
grep VQSRTrancheSNP /mnt/internal_2/priscilla/HSparents/variants/merged.hc.raw.vsqr.99.9.vcf > /mnt/internal_2/priscilla/HSparents/variants/merged.hc.raw.vsqr.tranch99.9.removedSNPs

for chr in 2L 2R 3L 3R X ; do
    echo $chr
    grep $chr /mnt/internal_2/priscilla/HSparents/variants/merged.hc.raw.vsqr.tranch99.9.removedSNPs | wc -l
done


#####
# END MORE LENIENT FILTERING
#####


#####
# CLEAN UP STRINGENT FILTERED VCF FOR HARP/RABBIT
#####



### remove repetitive regions from more stringent filtering and keep header  (wc = 3869453 )	
bedtools intersect -v \
	-header \
	-a /mnt/internal_2/priscilla/HSparents/variants/merged.hc.filter.vsqr.recode.vcf \
	-b /mnt/internal_2/priscilla/HSparents/variants/rm.edit.bed > \
	/mnt/internal_2/priscilla/HSparents/variants/merged.hc.filter.vsqr.recode.noREP.vcf

#annotate variants to sort out mnps etc
java -jar /usr/local/bin/GenomeAnalysisTK.jar \
   -R /mnt/internal_2/priscilla/ref_files/all_dmel.fasta \
   -T VariantAnnotator \
   -V /mnt/internal_2/priscilla/HSparents/variants/merged.hc.filter.vsqr.recode.noREP.vcf \
   -A VariantType \
   -o /mnt/internal_2/priscilla/HSparents/variants/merged.hc.filter.vsqr.recode.noREP.annotated.vcf


#remove MULTIALLELIC_SNPs (wc = 3778198)

grep -v MULTIALLELIC_SNP /mnt/internal_2/priscilla/HSparents/variants/merged.hc.filter.vsqr.recode.noREP.annotated.vcf > /mnt/internal_2/priscilla/HSparents/variants/merged.hc.filter.vsqr.recode.noREP.noMNP.annotated.vcf

#only work with the main chromosome arms (wc=3770760)

grep -e "##" -e "#CHROM" -e "2L" -e "2R" -e "3L" -e "3R" -e "X" /mnt/internal_2/priscilla/HSparents/variants/merged.hc.filter.vsqr.recode.noREP.noMNP.annotated.vcf > /mnt/internal_2/priscilla/HSparents/variants/merged.hc.filter.vsqr.recode.noREP.noMNP.annotated.main.vcf

#make a bed file of annotated indels +/- 100 bp (wc of indels = 799649 )
  
grep -e INSERTION -e DELETION -e MULTIALLELIC_COMPLEX -e MULTIALLELIC_MIXED /mnt/internal_2/priscilla/HSparents/variants/merged.hc.filter.vsqr.recode.noREP.noMNP.annotated.main.vcf | \
    cut -f1-2 | awk '{print $1"\t"$2-100"\t"$2+100}' > \
    /mnt/internal_2/priscilla/HSparents/variants/merged.hc.indel100.bed

#remove indels +/- 100 bp (wc = 1126916)
bedtools intersect -v \
	-header \
	-a /mnt/internal_2/priscilla/HSparents/variants/merged.hc.filter.vsqr.recode.noREP.noMNP.annotated.main.vcf \
	-b /mnt/internal_2/priscilla/HSparents/variants/merged.hc.indel100.bed > \
	/mnt/internal_2/priscilla/HSparents/variants/merged.hc.filter.vsqr.recode.noREP.noMNP.noINDEL.vcf

#try indels +/- 50bp  (wc= 1740783, saves extra half million snps)
  
#grep -e INSERTION -e DELETION -e MULTIALLELIC_COMPLEX -e MULTIALLELIC_MIXED /mnt/internal_2/priscilla/HSparents/variants/merged.hc.filter.vsqr.recode.noREP.noMNP.annotated.main.vcf | \
#    cut -f1-2 | awk '{print $1"\t"$2-50"\t"$2+50}' > \
#		    /mnt/internal_2/priscilla/HSparents/variants/merged.hc.indel50.bed

#remove indels +/- 50 bp (wc = )
#bedtools intersect -v \
#	-header \
#	-a /mnt/internal_2/priscilla/HSparents/variants/merged.hc.filter.vsqr.recode.noREP.noMNP.annotated.main.vcf \
#	-b /mnt/internal_2/priscilla/HSparents/variants/merged.hc.indel50.bed > \
#	/mnt/internal_2/priscilla/HSparents/variants/merged.hc.filter.vsqr.recode.noREP.noMNP.noINDEL50.vcf

#remove snps with median depth <10: (alan says nevermind)

#bedtools intersect -v \
#	-header \
#	-a /mnt/internal_2/priscilla/HSparents/variants/merged.hc.filter.vsqr.recode.noREP.noMNP.noINDEL.vcf \
#	-b /mnt/internal_2/priscilla/HSparents/variants/merged.hc.mediandepthbelow10.bed > \
#	/mnt/internal_2/priscilla/HSparents/variants/merged.hc.filter.vsqr.recode.noREP.noMNP.noINDEL.medDepth10.vcf


##########
# FILE CLEANUP AND SPLITTING
##########



#clean up vcf file to just have genotypes
/mnt/internal_2/priscilla/HSparents/Programs/vcflib/bin/vcfkeepgeno /mnt/internal_2/priscilla/HSparents/variants/merged.hc.filter.vsqr.recode.noREP.noMNP.noINDEL.vcf GT > /mnt/internal_2/priscilla/HSparents/variants/merged.hc.clean.vcf

#get rid of ugly info too
/mnt/internal_2/priscilla/HSparents/Programs/vcflib/bin/vcfkeepinfo /mnt/internal_2/priscilla/HSparents/variants/merged.hc.clean.vcf DP > /mnt/internal_2/priscilla/HSparents/variants/merged.hc.clean2.vcf


#need to make a SNP ID. cory says it is easiest to just parse out header and body, edit body, then rejoin

grep "#" /mnt/internal_2/priscilla/HSparents/variants/merged.hc.clean2.vcf > /mnt/internal_2/priscilla/HSparents/variants/merged.hc.clean2.header.vcf

grep -v "#" /mnt/internal_2/priscilla/HSparents/variants/merged.hc.clean2.vcf > /mnt/internal_2/priscilla/HSparents/variants/merged.hc.clean2.body.vcf

awk 'BEGIN { OFS = "\t"; } {$3 = $1"_"$2"_SNP"; print}'  /mnt/internal_2/priscilla/HSparents/variants/merged.hc.clean2.body.vcf >  /mnt/internal_2/priscilla/HSparents/variants/merged.hc.clean2.body.snpid.vcf

cat /mnt/internal_2/priscilla/HSparents/variants/merged.hc.clean2.header.vcf /mnt/internal_2/priscilla/HSparents/variants/merged.hc.clean2.body.snpid.vcf > /mnt/internal_2/priscilla/HSparents/variants/merged.hc.clean3.vcf

#make lists of A and B lines

#pull out the p### for lines A and B

awk '{if($3=="A") print $1 > "/mnt/internal_2/priscilla/HSparents/variants/A.txt"; else print $1 > "/mnt/internal_2/priscilla/HSparents/variants/B.txt"}' /mnt/pricey_2/priscilla/hybrid/etc/line_number_swarm.txt

#use this to subset into separate founders vcfs

bcftools view -S /mnt/internal_2/priscilla/HSparents/variants/A.txt /mnt/internal_2/priscilla/HSparents/variants/merged.hc.clean3.vcf > /mnt/pricey_2/priscilla/hybrid/seqs/merged.hc.clean.A.vcf

bcftools view -S /mnt/internal_2/priscilla/HSparents/variants/B.txt /mnt/internal_2/priscilla/HSparents/variants/merged.hc.clean3.vcf > /mnt/pricey_2/priscilla/hybrid/seqs/merged.hc.clean.B.vcf

# reheader files with p001-p034
bcftools reheader -s /mnt/pricey_2/priscilla/hybrid/etc/34lines.txt /mnt/pricey_2/priscilla/hybrid/seqs/merged.hc.clean.A.vcf > /mnt/pricey_2/priscilla/hybrid/seqs/hs.hc.A.vcf

bcftools reheader -s /mnt/pricey_2/priscilla/hybrid/etc/34lines.txt /mnt/pricey_2/priscilla/hybrid/seqs/merged.hc.clean.B.vcf > /mnt/pricey_2/priscilla/hybrid/seqs/hs.hc.B.vcf

#make files with old and new name keys

paste <(bcftools query -l /mnt/pricey_2/priscilla/hybrid/seqs/merged.hc.clean.A.vcf) <(cat /mnt/pricey_2/priscilla/hybrid/etc/34lines.txt) > /mnt/pricey_2/priscilla/hybrid/seqs/parent_ids.hc.A.txt
paste <(bcftools query -l /mnt/pricey_2/priscilla/hybrid/seqs/merged.hc.clean.B.vcf) <(cat /mnt/pricey_2/priscilla/hybrid/etc/34lines.txt) > /mnt/pricey_2/priscilla/hybrid/seqs/parent_ids.hc.B.txt

#move haplotype caller files to new directory to avoid confusion

mv /mnt/pricey_2/priscilla/hybrid/seqs/*hc* /mnt/pricey_2/priscilla/hybrid/hc

#prepare all files need for reconstruction

for swarm in A B; do
    
    #remove ## lines from swarm vcfs in order to remake priors files
    grep -v '##' /mnt/pricey_2/priscilla/hybrid/hc/hs.hc.$swarm.vcf > /mnt/pricey_2/priscilla/hybrid/hc/hs.hc.$swarm.noheader.vcf

done

# make priors files

/mnt/pricey_2/priscilla/hybrid/scripts/make_harp_input4_byswarm.sh

##define intervalList function

intervalList () {
fn=${1}
chr=$( head -n1 ${fn} | cut -f1 -d',' )

tail -n +2 ${fn} | cut -f1 -d',' | sed 's/\([0-9]*\)/\1-\1/g' | sed "s/^/${chr}:/g"  > ${fn}.intervals
	
}
export -f intervalList



#deal with priors

for swarm in A B; do
    
    for chr in 2L 2R 3L 3R X; do
	
	mkdir -p /mnt/pricey_2/priscilla/hybrid/hc/"$chr"
	mkdir -p /mnt/pricey_2/priscilla/hybrid/hc/"$chr"/"$swarm"

	#move prior to correct folder
	mv /mnt/pricey_2/priscilla/hybrid/hc/prior_"$swarm"_"$chr".noRep.noIndel.noHet.csv /mnt/pricey_2/priscilla/hybrid/hc/"$chr"/"$swarm"

	#index the snp table
	/usr/bin/index_snp_table \
	/mnt/pricey_2/priscilla/hybrid/hc/$chr/$swarm/prior_"$swarm"_"$chr".noRep.noIndel.noHet.csv \
	10000

	#make interval list

	intervalList /mnt/pricey_2/priscilla/hybrid/hc/$chr/$swarm/prior_"$swarm"_"$chr".noRep.noIndel.noHet.csv
	done
done


#make founder files

for swarm in A B; do
    
    for chr in 2L 2R 3L 3R X; do

	cd /mnt/pricey_2/priscilla/hybrid/hc/$chr/$swarm

	## pull out single chromosome vcfs
	grep  -e "#" -e $chr /mnt/pricey_2/priscilla/hybrid/hc/hs.hc."$swarm".vcf > /mnt/pricey_2/priscilla/hybrid/hc/$chr/$swarm/hs.hc.$chr.$swarm.vcf

	## make founder vcf
	/mnt/pricey_2/priscilla/hybrid/scripts/prepareFounderVCF.sh /mnt/pricey_2/priscilla/hybrid/hc/$chr/$swarm/hs.hc.$chr.$swarm.vcf "$chr"

	## build founder matrix
	Rscript /mnt/pricey_2/priscilla/hybrid/scripts/buildFounderMatrix_generalized.Rscript /mnt/pricey_2/priscilla/hybrid/hc/$chr/$swarm/hs.hc.$chr.$swarm.vcf.replaced /mnt/pricey_2/priscilla/hybrid/etc/recombination.bed "$chr" /mnt/pricey_2/priscilla/hybrid/scripts

	## make rabbit.index
	awk -v chromosome="$chr" 'chromosome ~ $1 {print $2}' /mnt/pricey_2/priscilla/hybrid/hc/$chr/$swarm/hs.hc.$chr.$swarm.vcf | awk 'BEGIN{OFS="\t";} {print NR,$0}' > /mnt/pricey_2/priscilla/hybrid/hc/$chr/$swarm/RABBIT.index

    done
done


	
#get everythign to rivanna
	

rsync -avzhe ssh /mnt/pricey_2/priscilla/hybrid/hc pae3g@rivanna1.hpc.virginia.edu:/scratch/pae3g/hybrid



#make file that includes filtered indels for reconstruction

#######
#INDEL VQSR to make an all polymorphisms file for genome reconstruction
#######

#get dgrp indels (wc=462423)
head -n 7 /mnt/internal_2/priscilla/HSparents/variants/dgrp2.vcf > /mnt/internal_2/priscilla/HSparents/variants/dgrp2_indel.vcf     
grep -E "INS|DEL" /mnt/internal_2/priscilla/HSparents/variants/dgrp2.vcf >> /mnt/internal_2/priscilla/HSparents/variants/dgrp2_indel.vcf

#downsample indels (wc=23456)
/mnt/internal_2/priscilla/vcflib/bin/vcfrandomsample -r 0.05 -p 420  /mnt/internal_2/priscilla/HSparents/variants/dgrp2_indel.vcf > /mnt/internal_2/priscilla/HSparents/variants/dgrp2_indel_randomSubset.vcf

#intersect dgrp indels with mine to make indel truthiness set (wc=

bedtools intersect -header \
	 -a /mnt/internal_2/priscilla/HSparents/variants/dgrp2_indel_randomSubset.vcf \
	 -b /mnt/internal_2/priscilla/HSparents/variants/merged.hc.raw.vcf > \
	 /mnt/internal_2/priscilla/HSparents/variants/indel_dgrp2_randomSubset.truthiness.vcf


#run recalibration with this truthiness file on file with less stringent recal applied

java -jar /usr/local/bin/GenomeAnalysisTK.jar \
      -T VariantRecalibrator \
      -R /mnt/internal_2/priscilla/ref_files/all_dmel.fasta \
      -input /mnt/internal_2/priscilla/HSparents/variants/merged.hc.raw.vsqr.99.9.vcf   \
      -resource:dgrp,known=false,training=true,truth=true,prior=15.0 /mnt/internal_2/priscilla/HSparents/variants/indel_dgrp2_randomSubset.truthiness.vcf \
      -an DP \
      -an QD \
      -an FS \
      -an SOR \
      -an MQ \
      -an MQRankSum \
      -an ReadPosRankSum \
      -mode INDEL \
      -nt 12 \
      -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
      -recalFile /mnt/internal_2/priscilla/HSparents/variants/hc.recalibrate_INDEL_merged.recal \
      -tranchesFile /mnt/internal_2/priscilla/HSparents/variants/hc.recalibrate_INDEL.tranches \
      -rscriptFile /mnt/internal_2/priscilla/HSparents/variants/hc.recalibrate_INDEL_plots_merged.R 

#apply recalibration at a more lenient threshold to increase snp count 
java -jar /usr/local/bin/GenomeAnalysisTK.jar \
      -T ApplyRecalibration \
      -R /mnt/internal_2/priscilla/ref_files/all_dmel.fasta \
      -input /mnt/internal_2/priscilla/HSparents/variants/merged.hc.raw.vsqr.99.9.vcf \
      -mode INDEL \
      --ts_filter_level 99.9 \
      -recalFile /mnt/internal_2/priscilla/HSparents/variants/hc.recalibrate_INDEL_merged.recal \
      -tranchesFile /mnt/internal_2/priscilla/HSparents/variants/hc.recalibrate_INDEL.tranches \
      -o /mnt/internal_2/priscilla/HSparents/variants/merged.hc.raw.vsqr.snpindel.99.9.vcf  


#clean up with filter with more lenient 99.9% tranch (keeps 4400410 variants)
vcftools \
--vcf /mnt/internal_2/priscilla/HSparents/variants/merged.hc.raw.vsqr.snpindel.99.9.vcf  \
--remove-filtered-all \
--recode \
--recode-INFO-all \
--out /mnt/internal_2/priscilla/HSparents/variants/merged.hc.filter.vsqr.snpindel.99.9

### remove repetitive regions from less stringent filtering and keep header  (~4.1 million snps after)	
bedtools intersect -v \
	-header \
	-a /mnt/internal_2/priscilla/HSparents/variants/merged.hc.filter.vsqr.snpindel.99.9.recode.vcf \
	-b /mnt/internal_2/priscilla/HSparents/variants/rm.edit.bed > \
	/mnt/internal_2/priscilla/HSparents/variants/merged.hc.filter.vsqr.snpindel.99.9.recode.noREP.vcf

#annotate variants to sort out mnps etc
java -jar /usr/local/bin/GenomeAnalysisTK.jar \
   -R /mnt/internal_2/priscilla/ref_files/all_dmel.fasta \
   -T VariantAnnotator \
   -V /mnt/internal_2/priscilla/HSparents/variants/merged.hc.filter.vsqr.snpindel.99.9.recode.noREP.vcf \
   -A VariantType \
   -o /mnt/internal_2/priscilla/HSparents/variants/merged.hc.filter.vsqr.snpindel.99.9.recode.noREP.annotated.vcf

#remove multiallelic sites from less stringent filtering (wc = 3.7 million)

grep -v MULTIALLELIC /mnt/internal_2/priscilla/HSparents/variants/merged.hc.filter.vsqr.snpindel.99.9.recode.noREP.annotated.vcf > /mnt/internal_2/priscilla/HSparents/variants/merged.hc.filter.vsqr.snpindel.99.9.recode.noREP.biallelic.vcf

#####
#CLEAN UP TO MATCH OTHER VCFS--NUMERIC PARENTS AND GENOTYPES ONLY, SPLIT BY SWARM AND CHR
#####


#clean up and split by swarm and chromosome:

/mnt/internal_2/priscilla/HSparents/Programs/vcflib/bin/vcfkeepgeno /mnt/internal_2/priscilla/HSparents/variants/merged.hc.filter.vsqr.snpindel.99.9.recode.noREP.biallelic.vcf GT > /mnt/internal_2/priscilla/HSparents/variants/merged.hc.snp.indel.clean.vcf

#get rid of ugly info too
/mnt/internal_2/priscilla/HSparents/Programs/vcflib/bin/vcfkeepinfo /mnt/internal_2/priscilla/HSparents/variants/merged.hc.snp.indel.clean.vcf VariantType > /mnt/internal_2/priscilla/HSparents/variants/merged.hc.snp.indel.clean2.vcf



#need to make a VARIANT ID. cory says it is easiest to just parse out header and body, edit body, then rejoin

grep "#" /mnt/internal_2/priscilla/HSparents/variants/merged.hc.snp.indel.clean2.vcf > /mnt/internal_2/priscilla/HSparents/variants/merged.hc.snp.indel.clean2.header.vcf

grep -v "#" /mnt/internal_2/priscilla/HSparents/variants/merged.hc.snp.indel.clean2.vcf > /mnt/internal_2/priscilla/HSparents/variants/merged.hc.snp.indel.clean2.body.vcf


awk 'BEGIN { OFS = "\t"; } {$3 = $1"_"$2"_VAR"; print}'  /mnt/internal_2/priscilla/HSparents/variants/merged.hc.snp.indel.clean2.body.vcf >  /mnt/internal_2/priscilla/HSparents/variants/merged.hc.snp.indel.clean2.body.snpid.vcf

cat /mnt/internal_2/priscilla/HSparents/variants/merged.hc.snp.indel.clean2.header.vcf /mnt/internal_2/priscilla/HSparents/variants/merged.hc.snp.indel.clean2.body.snpid.vcf > /mnt/internal_2/priscilla/HSparents/variants/merged.hc.snp.indel.clean3.vcf

#make lists of A and B lines

#use this to subset into separate founders vcfs

bcftools view -S /mnt/internal_2/priscilla/HSparents/variants/A.txt /mnt/internal_2/priscilla/HSparents/variants/merged.hc.snp.indel.clean3.vcf > /mnt/pricey_2/priscilla/hybrid/hc/merged.hc.snp.indel.clean.A.vcf

bcftools view -S /mnt/internal_2/priscilla/HSparents/variants/B.txt /mnt/internal_2/priscilla/HSparents/variants/merged.hc.snp.indel.clean3.vcf > /mnt/pricey_2/priscilla/hybrid/hc/merged.hc.snp.indel.clean.B.vcf

# reheader files with p001-p034
bcftools reheader -s /mnt/pricey_2/priscilla/hybrid/etc/34lines.txt /mnt/pricey_2/priscilla/hybrid/hc/merged.hc.snp.indel.clean.A.vcf > /mnt/pricey_2/priscilla/hybrid/hc/hs.hc.snp.indel.A.vcf

bcftools reheader -s /mnt/pricey_2/priscilla/hybrid/etc/34lines.txt /mnt/pricey_2/priscilla/hybrid/hc/merged.hc.snp.indel.clean.B.vcf > /mnt/pricey_2/priscilla/hybrid/hc/hs.hc.snp.indel.B.vcf

for swarm in A B; do
    
    for chr in 2L 2R 3L 3R X; do

	cd /mnt/pricey_2/priscilla/hybrid/hc/$chr/$swarm

	## pull out single chromosome vcfs
	grep  -e "#" -e $chr /mnt/pricey_2/priscilla/hybrid/hc/hs.hc.snp.indel."$swarm".vcf > /mnt/pricey_2/priscilla/hybrid/hc/$chr/$swarm/hs.hc.snp.indel.$chr.$swarm.vcf
 
    done
done


#test out vcf for reconstruction

/mnt/pricey_2/priscilla/hybrid/scripts/tabixIndexRealData_generalized.py "PAE031_H9" "2L" hs.hc.snp.indel.2L.B.vcf /mnt/pricey_2/priscilla/hybrid//etc/htslib "round2.estimate.paths" 


#get reconstructed vcfs off of rivanna

nohup rsync -avzh pae3g@rivanna1.hpc.virginia.edu:/scratch/pae3g/hybrid/hc/2L/A/archive_10founders10Ksnps /mnt/pricey_2/priscilla/hybrid/hc/2L/A &

nohup rsync -avzh pae3g@rivanna1.hpc.virginia.edu:/scratch/pae3g/hybrid/hc/2L/A/VCFs_10founders /mnt/pricey_2/priscilla/hybrid/hc/2L/A &



nohup rsync -avzh pae3g@rivanna1.hpc.virginia.edu:/scratch/pae3g/hybrid/hc/2L/B/archive /mnt/pricey_2/priscilla/hybrid/hc/2L/B &

nohup rsync -avzh pae3g@rivanna1.hpc.virginia.edu:/scratch/pae3g/hybrid/hc/2L/B/VCFs /mnt/pricey_2/priscilla/hybrid/hc/2L/B &




########
#MAKE ANOTHER FILE THAT DOESN'T HAVE INDELS
########

#remove INDEL sites from less stringent filtering (wc = 3.1 million )

grep -v INSERTION /mnt/internal_2/priscilla/HSparents/variants/merged.hc.filter.vsqr.snpindel.99.9.recode.noREP.biallelic.vcf | grep -v DELETION > /mnt/internal_2/priscilla/HSparents/variants/merged.hc.filter.vsqr.snponly.99.9.recode.noREP.biallelic.vcf


#clean up and split by swarm and chromosome:

/mnt/internal_2/priscilla/HSparents/Programs/vcflib/bin/vcfkeepgeno /mnt/internal_2/priscilla/HSparents/variants/merged.hc.filter.vsqr.snponly.99.9.recode.noREP.biallelic.vcf GT > /mnt/internal_2/priscilla/HSparents/variants/merged.hc.snp.99.9.clean.vcf

#get rid of ugly info too
/mnt/internal_2/priscilla/HSparents/Programs/vcflib/bin/vcfkeepinfo /mnt/internal_2/priscilla/HSparents/variants/merged.hc.snp.99.9.clean.vcf VariantType > /mnt/internal_2/priscilla/HSparents/variants/merged.hc.snp.99.9.clean2.vcf



#need to make a VARIANT ID. cory says it is easiest to just parse out header and body, edit body, then rejoin

grep "#" /mnt/internal_2/priscilla/HSparents/variants/merged.hc.snp.99.9.clean2.vcf > /mnt/internal_2/priscilla/HSparents/variants/merged.hc.snp.99.9.clean2.header.vcf

grep -v "#" /mnt/internal_2/priscilla/HSparents/variants/merged.hc.snp.99.9.clean2.vcf > /mnt/internal_2/priscilla/HSparents/variants/merged.hc.snp.99.9.clean2.body.vcf


awk 'BEGIN { OFS = "\t"; } {$3 = $1"_"$2"_SNP"; print}'  /mnt/internal_2/priscilla/HSparents/variants/merged.hc.snp.99.9.clean2.body.vcf >  /mnt/internal_2/priscilla/HSparents/variants/merged.hc.snp.99.9.clean2.body.snpid.vcf

cat /mnt/internal_2/priscilla/HSparents/variants/merged.hc.snp.99.9.clean2.header.vcf /mnt/internal_2/priscilla/HSparents/variants/merged.hc.snp.99.9.clean2.body.snpid.vcf > /mnt/internal_2/priscilla/HSparents/variants/merged.hc.snp.99.9.clean3.vcf

#make lists of A and B lines

#use this to subset into separate founders vcfs

bcftools view -S /mnt/internal_2/priscilla/HSparents/variants/A.txt /mnt/internal_2/priscilla/HSparents/variants/merged.hc.snp.99.9.clean3.vcf > /mnt/pricey_2/priscilla/hybrid/hc/merged.hc.snp.99.9.clean.A.vcf

bcftools view -S /mnt/internal_2/priscilla/HSparents/variants/B.txt /mnt/internal_2/priscilla/HSparents/variants/merged.hc.snp.99.9.clean3.vcf > /mnt/pricey_2/priscilla/hybrid/hc/merged.hc.snp.99.9.clean.B.vcf

# reheader files with p001-p034
bcftools reheader -s /mnt/pricey_2/priscilla/hybrid/etc/34lines.txt /mnt/pricey_2/priscilla/hybrid/hc/merged.hc.snp.99.9.clean.A.vcf > /mnt/pricey_2/priscilla/hybrid/hc/hs.hc.snp.99.9.A.vcf

bcftools reheader -s /mnt/pricey_2/priscilla/hybrid/etc/34lines.txt /mnt/pricey_2/priscilla/hybrid/hc/merged.hc.snp.99.9.clean.B.vcf > /mnt/pricey_2/priscilla/hybrid/hc/hs.hc.snp.99.9.B.vcf

for swarm in A B; do
    
    for chr in 2L 2R 3L 3R X; do

	cd /mnt/pricey_2/priscilla/hybrid/hc/$chr/$swarm

	## pull out single chromosome vcfs
	grep -v "##" /mnt/pricey_2/priscilla/hybrid/hc/hs.hc.snp.99.9."$swarm".vcf | grep -e "#" -e $chr  > /mnt/pricey_2/priscilla/hybrid/hc/$chr/$swarm/hs.hc.snp.99.9.$chr.$swarm.vcf
 
    done
done

