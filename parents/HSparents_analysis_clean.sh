
#Hybrid swarm parents HiSeqX10 lane

#names of all files are in in "filenames_HKCHFALXX_sorted.txt"
#all sample information is in "HKCHFALXX_all_file_info.txt"

#make a text file that has all of the file stems to be used for mapping parallel
cut -f5 filenames_HKCHFALXX_sorted.txt | awk  '{print $2}' FS='s7_[1-2]_' | uniq > samples.txt

#run script that includes parallel command to do all mapping
/mnt/internal_2/priscilla/HSparents/scripts/pear_bwa.sh

#map all SRA reads
/mnt/internal_2/priscilla/HSparents/scripts/map_sra2.sh

#make a list of all bam files

ls /mnt/internal_2/priscilla/HSparents/mapped/*.sort.dedup.bam > /mnt/internal_2/priscilla/HSparents/scripts/bams.list

#use unified genotyper to make an initial VCF of new sequence data

nohup java -jar /usr/local/bin/GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-nt 4 \
-nct 6 \
-hets 0.01 \
-R /mnt/internal_2/priscilla/ref_files/all_dmel.fasta \
-I /mnt/internal_2/priscilla/HSparents/scripts/bams.list  \
-o /mnt/internal_2/priscilla/HSparents/variants/hsparents.ug.raw.vcf &

#this file is a random subset of known DGRP SNPs to serve as "truth" SNPs (available on dataDryad)
inbred.ug.dgrp2_snp.randomSubset.noREP.noINDEL.vcf.truthiness.vcf 

#this file has tranch cutoffs for SNP verification
recalibrate_SNP.tranches 

#this bed file is used to mask all repeats +/- 100 bp
repMasker.dgrp2_indel_100bp.sort.bed 

 ### run VariantRecalibrator
nohup java -jar /usr/local/bin/GenomeAnalysisTK.jar \
     -T VariantRecalibrator \
-R /mnt/internal_2/priscilla/ref_files/all_dmel.fasta \
-input /mnt/internal_2/priscilla/HSparents/variants/hsparents.ug.raw.vcf \
-resource:dgrp,known=false,training=true,truth=true,prior=15.0 /mnt/internal_2/priscilla/HSparents/variants/inbred.ug.dgrp2_snp.randomSubset.noREP.noINDEL.vcf.truthiness.vcf \
-an DP \
-an QD \
-an FS \
-an SOR \
-an MQ \
-an MQRankSum \
-an ReadPosRankSum \
-mode SNP \
-nt 10 \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
-recalFile /mnt/internal_2/priscilla/HSparents/variants/recalibrate_SNP.recal \
-tranchesFile /mnt/internal_2/priscilla/HSparents/variants/recalibrate_SNP.tranches \
-rscriptFile /mnt/internal_2/priscilla/HSparents/variants/recalibrate_SNP_plots.R &



### run ApplyRecalibration
nohup java -jar /usr/local/bin/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R /mnt/internal_2/priscilla/ref_files/all_dmel.fasta \
-input /mnt/internal_2/priscilla/HSparents/variants/hsparents.ug.raw.vcf \
-mode SNP \
--ts_filter_level 99.0 \
-recalFile /mnt/internal_2/priscilla/HSparents/variants/recalibrate_SNP.recal \
-tranchesFile /mnt/internal_2/priscilla/HSparents/variants/recalibrate_SNP.tranches \
-o /mnt/internal_2/priscilla/HSparents/variants/hsparents.ug.raw.vsqr.vcf &

### filter sites based on PASS field
vcftools \
--vcf /mnt/internal_2/priscilla/HSparents/variants/hsparents.ug.raw.vsqr.vcf \
--remove-filtered-all \
--recode \
--recode-INFO-all \
--out /mnt/internal_2/priscilla/HSparents/variants/hsparents.ug.filter.vsqr

### remove repetitive regions & filter around indels & convert to BCF
bedtools intersect -sorted -v -header \
-a /mnt/internal_2/priscilla/HSparents/variants/hsparents.ug.filter.vsqr.recode.vcf \
-b /mnt/internal_2/priscilla/HSparents/variants/repMasker.dgrp2_indel_100bp.sort.bed | \
tee /mnt/internal_2/priscilla/HSparents/variants/hsparents.ug.filter.vsqr.recode.noREP.noINDEL.vcf | \
bgzip -c > \
/mnt/internal_2/priscilla/HSparents/variants/hsparents.filter.vsqr.recode.noREP.noINDEL.vcf.gz

#get list of sample names in each vcf
bcftools query -l inbred.ug.filter.vsqr.recode.noREP.noINDEL.vcf.gz >inbred_samples.txt
bcftools query -l hsparents.ug.filter.vsqr.recode.noREP.noINDEL.vcf > hsparents_samples.txt

cp inbred_samples.txt inbred_samples_edit.txt

sed -i 's/Ithica_//g' inbred_samples_edit.txt
sed -i 's/\-/\_/g' inbred_samples_edit.txt
sed -i 's/(//g' inbred_samples_edit.txt
sed -i 's/)//g' inbred_samples_edit.txt
sed -i 's/Birmingham_AL_//g' inbred_samples_edit.txt
sed -i 's/Selba_AL_//g' inbred_samples_edit.txt
sed -i 's/TampaBay_FL_//g' inbred_samples_edit.txt
sed -i 's/Thomasville_GA_//g' inbred_samples_edit.txt
sed -i 's/BullocksHarbor_BerryIslands__//g' inbred_samples_edit.txt
sed -i 's/Cockburn_SanSalvador_4//g' inbred_samples_edit.txt
sed -i 's/Freeport_GrandBahamasWest_//g' inbred_samples_edit.txt
sed -i 's/GeorgeTown_Exumas_//g' inbred_samples_edit.txt
sed -i 's/Mayaguana_Mayaguana_//g' inbred_samples_edit.txt
sed -i 's/Meridian_MS_//g' inbred_samples_edit.txt

Paste inbred_samples.txt inbred_samples_edit.txt > inbred_samples_both_names.txt

tail -n +2 inbred_samples_both_names.txt | sort -k 2 > inbred_samples_sorted.txt

Sort hsparents_samples.txt > hsparents_samples_sorted.txt

join hsparents_samples_sorted.txt inbred_samples_sorted.txt -2 2 > samples_used_both_names.txt



#Going to remake ‘truth’ VCF to include missing DGRP files. First move all relevant bam files to a new folder (inbred_mapped) by looping through a file list and grepping all.

awk '{print $2}' samples_used_both_names.txt > bam_files_to_move.txt

while read pattern; do 
    echo cp /mnt/icy_1/inbredLines/01_mapped/*${pattern}*  /mnt/internal_2/priscilla/HSparents/inbred_mapped/ done < bam_files_to_move.txt

while read pattern; do cp /mnt/icy_1/inbredLines/01_mapped/*${pattern}* /mnt/internal_2/priscilla/HSparents/inbred_mapped; done < bam_files_to_move.txt 
#that didn’t get everything due to different punctuation, so also do this:

while read pattern; do cp /mnt/icy_1/inbredLines/01_mapped/*${pattern}* /mnt/internal_2/priscilla/HSparents/inbred_mapped; done < hsparents_samples_sorted.txt 


#Download DGRP data

nohup sh -c 'while read srx strain sra; do /mnt/internal_2/priscilla/HSparents/Programs/sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump --outdir /mnt/internal_2/priscilla/HSparents/SRAdata/$strain --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-files --clip $sra; done < sra_dgrp.txt' &


#redo DGRP download because bwa doesn’t like the read ids that were appended

nohup sh -c 'while read srx strain sra; do 
/mnt/internal_2/priscilla/HSparents/Programs/sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump --outdir /mnt/internal_2/priscilla/HSparents/SRAdata/$strain --gzip --skip-technical  --read-filter pass --dumpbase --split-files --clip $sra; done < sra_dgrp.txt' &

#map DGRP SRA data
nohup /mnt/internal_2/priscilla/HSparents/scripts/map_sra.sh &
#PID 26318

#DGRP28124 is actually single end data, so will need to redo download with correct parameters


nohup /mnt/internal_2/priscilla/HSparents/Programs/sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump --outdir /mnt/internal_2/priscilla/HSparents/SRAdata/DGRP_28124 --gzip --skip-technical  --read-filter pass --dumpbase --clip SRR933567 &

./map_sra2.sh #new version to deal with single end data names and freshly dowloaded data




#Alyssa moved merged bam files for 6 of the 12BME10 lines from Rivanna to /mnt/internal_2/tempInbred
cp /mnt/internal_2/tempInbred/* /mnt/internal_2/priscilla/HSparents/inbred_mapped

ls *merge.bam > ../scripts/bam1.txt
ls *up.bam > ../scripts/bam2.txt
cat bam1.txt bam2.txt > inbred_bams.list



nohup java -jar /usr/local/bin/GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-nt 4 \
-nct 6 \
-hets 0.01 \
-R /mnt/internal_2/priscilla/ref_files/all_dmel.fasta \
-I /mnt/internal_2/priscilla/HSparents/scripts/inbred_bams.list  \
-o /mnt/internal_2/priscilla/HSparents/variants/inbred_68.ug.raw.vcf &
#PID 24774

 ### run VariantRecalibrator
nohup java -jar /usr/local/bin/GenomeAnalysisTK.jar \
     -T VariantRecalibrator \
-R /mnt/internal_2/priscilla/ref_files/all_dmel.fasta \
-input /mnt/internal_2/priscilla/HSparents/variants/inbred_68.ug.raw.vcf \
-resource:dgrp,known=false,training=true,truth=true,prior=15.0 /mnt/internal_2/priscilla/HSparents/variants/inbred.ug.dgrp2_snp.randomSubset.noREP.noINDEL.vcf.truthiness.vcf \
-an DP \
-an QD \
-an FS \
-an SOR \
-an MQ \
-an MQRankSum \
-an ReadPosRankSum \
-mode SNP \
-nt 10 \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
-recalFile /mnt/internal_2/priscilla/HSparents/variants/recalibrate_SNP_inbred_68.recal \
-tranchesFile /mnt/internal_2/priscilla/HSparents/variants/recalibrate_SNP.tranches \
-rscriptFile /mnt/internal_2/priscilla/HSparents/variants/recalibrate_SNP_plots_inbred_68.R &

#PID 25113

nohup java -jar /usr/local/bin/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R /mnt/internal_2/priscilla/ref_files/all_dmel.fasta \
-input /mnt/internal_2/priscilla/HSparents/variants/inbred_68.ug.raw.vcf \
-mode SNP \
--ts_filter_level 99.0 \
-recalFile /mnt/internal_2/priscilla/HSparents/variants/recalibrate_SNP_inbred_68.recal \
-tranchesFile /mnt/internal_2/priscilla/HSparents/variants/recalibrate_SNP.tranches \
-o /mnt/internal_2/priscilla/HSparents/variants/inbred_68.ug.raw.vsqr.vcf &

#PID  29964

### filter sites based on PASS field
vcftools \
--vcf /mnt/internal_2/priscilla/HSparents/variants/inbred_68.ug.raw.vsqr.vcf \
--remove-filtered-all \
--recode \
--recode-INFO-all \
--out /mnt/internal_2/priscilla/HSparents/variants/inbred_68.ug.filter.vsqr

### remove repetitive regions & filter around indels & convert to BCF
bedtools intersect -sorted -v -header \
-a /mnt/internal_2/priscilla/HSparents/variants/inbred_68.ug.filter.vsqr.recode.vcf \
-b /mnt/internal_2/priscilla/HSparents/variants/repMasker.dgrp2_indel_100bp.sort.bed | \
tee /mnt/internal_2/priscilla/HSparents/variants/inbred_68.ug.filter.vsqr.recode.noREP.noINDEL.vcf | \
bgzip -c > \
/mnt/internal_2/priscilla/HSparents/variants/inbred_68.ug.filter.vsqr.recode.noREP.noINDEL.vcf.gz

bcftools query -l inbred_68.ug.filter.vsqr.recode.noREP.noINDEL.vcf.gz >inbred_68_samples.txt

cp inbred_68_samples.txt inbred_68_samples_edit.txt

sed -i 's/Ithica_//g' inbred_68_samples_edit.txt
sed -i 's/\-/\_/g' inbred_68_samples_edit.txt
sed -i 's/(//g' inbred_68_samples_edit.txt
sed -i 's/)//g' inbred_68_samples_edit.txt
sed -i 's/Birmingham_AL_//g' inbred_68_samples_edit.txt
sed -i 's/Selba_AL_//g' inbred_68_samples_edit.txt
sed -i 's/TampaBay_FL_//g' inbred_68_samples_edit.txt
sed -i 's/Thomasville_GA_//g' inbred_68_samples_edit.txt
sed -i 's/BullocksHarbor_BerryIslands_//g' inbred_68_samples_edit.txt
sed -i 's/Cockburn_SanSalvador_4//g' inbred_68_samples_edit.txt
sed -i 's/Freeport_GrandBahamasWest_//g' inbred_68_samples_edit.txt
sed -i 's/GeorgeTown_Exumas_//g' inbred_68_samples_edit.txt
sed -i 's/Mayaguana_Mayaguana_//g' inbred_68_samples_edit.txt
sed -i 's/Meridian_MS_//g' inbred_68_samples_edit.txt

#make copy of inbred_68 vcf file
cp inbred_68.ug.filter.vsqr.recode.noREP.noINDEL.vcf.gz inbred_68.ug.filter.vsqr.recode.noREP.noINDEL_copy.vcf.gz

#attempt to rename samples in bcftools

bcftools reheader -s inbred_68_samples_edit.txt inbred_68.ug.filter.vsqr.recode.noREP.noINDEL_copy.vcf.gz > inbred_68.ug.filter.vsqr.recode.noREP.noINDEL_renamed.vcf.gz

#unzip and keep a zipped copy
gunzip -k inbred_68.ug.filter.vsqr.recode.noREP.noINDEL_renamed.vcf.gz

#Split each vcf into separate vcfs to use gtcompare tool;


bcftools query -l inbred_68.ug.filter.vsqr.recode.noREP.noINDEL_renamed.vcf | parallel -j 8 'vcf-subset --exclude-ref -c {} inbred_68.ug.filter.vsqr.recode.noREP.noINDEL_renamed.vcf > {}.inbred.vcf'


bcftools query -l hsparents.ug.filter.vsqr.recode.noREP.noINDEL_renamed.vcf | parallel -j 8 'vcf-subset --exclude-ref -c {} hsparents.ug.filter.vsqr.recode.noREP.noINDEL_renamed.vcf > {}.hsparent.vcf'

for sample in `bcftools query -l inbred_68.ug.filter.vsqr.recode.noREP.noINDEL_renamed.vcf`
do 
vcf-subset --exclude-ref -c $sample inbred_68.ug.filter.vsqr.recode.noREP.noINDEL_renamed.vcf > ${sample}.inbred.vcf
done



/mnt/internal_2/priscilla/HSparents/Programs/vcflib/bin/vcfannotategenotypes inbred 12BME10_104.hsparent_copy.vcf 12BME10_104.inbred_copy.vcf > 12BME10_104_annotated.vcf


/mnt/internal_2/priscilla/HSparents/Programs/vcflib/bin/vcfgenotypecompare inbred 12BME10_104_annotated.vcf > 12BME10_104_gtcompare.vcf

/mnt/internal_2/priscilla/HSparents/Programs/vcflib/bin/vcf2tsv 12BME10_104_gtcompare.vcf > 12BME10_104.tsv

/mnt/internal_2/priscilla/HSparents/Programs/vcflib/bin/vcf2tsv -n NA -g 12BME10_104_gtcompare.vcf > 12BME10_104g.tsv


#In R:
this_other <- read.delim("12BME10_104g.tsv")
this_other[is.na(this_other)] <- 0 


gterror <- function(s, n) {
   with(s,
   (sum(eval(as.symbol(paste(n, ".genotypes.AA_AR", sep=""))))
    + sum(eval(as.symbol(paste(n, ".genotypes.AR_AA", sep=""))))
    + sum(eval(as.symbol(paste(n, ".genotypes.AA_RR", sep=""))))
    + sum(eval(as.symbol(paste(n, ".genotypes.RR_AA", sep=""))))
    + sum(eval(as.symbol(paste(n, ".genotypes.AR_RR", sep=""))))
    + sum(eval(as.symbol(paste(n, ".genotypes.RR_AR", sep="")))))
   / ((sum(eval(as.symbol(paste(n, ".genotypes.RR_RR", sep=""))))
       + sum(eval(as.symbol(paste(n, ".genotypes.AA_AA", sep=""))))
       + sum(eval(as.symbol(paste(n, ".genotypes.AR_AR", sep="")))))
      + (sum(eval(as.symbol(paste(n, ".genotypes.AA_AR", sep=""))))
         + sum(eval(as.symbol(paste(n, ".genotypes.AR_AA", sep=""))))
         + sum(eval(as.symbol(paste(n, ".genotypes.AA_RR", sep=""))))
         + sum(eval(as.symbol(paste(n, ".genotypes.RR_AA", sep=""))))
         + sum(eval(as.symbol(paste(n, ".genotypes.AR_RR", sep=""))))
         + sum(eval(as.symbol(paste(n, ".genotypes.RR_AR", sep="")))))))
}
gterror(this_other, "inbred") 

/mnt/internal_2/priscilla/HSparents/Programs/vcflib/bin/vcfannotategenotypes inbred 12BME10_108.hsparent.vcf 12BME10_108.inbred.vcf >  12BME10_108_annotated.vcf
/mnt/internal_2/priscilla/HSparents/Programs/vcflib/bin/vcfgenotypecompare inbred 12BME10_108_annotated.vcf > 12BME10_108_gtcompare.vcf
/mnt/internal_2/priscilla/HSparents/Programs/vcflib/bin/vcf2tsv -n NA -g 12BME10_108_gtcompare.vcf > 12BME10_108.tsv


#can you compare all genotypes at once?
nohup /mnt/internal_2/priscilla/HSparents/Programs/vcflib/bin/vcfannotategenotypes inbred hsparents.ug.filter.vsqr.recode.noREP.noINDEL.vcf inbred_68.ug.filter.vsqr.recode.noREP.noINDEL_renamed.vcf >  genotypes_annotated.vcf &

nohup /mnt/internal_2/priscilla/HSparents/Programs/vcflib/bin/vcfgenotypecompare inbred genotypes_annotated.vcf > genotypes_gtcompare.vcf &
#PID 30397

/mnt/internal_2/priscilla/HSparents/Programs/vcflib/bin/vcf2tsv -n NA -g genotypes_gtcompare.vcf > genotypes.tsv

nohup mnt/internal_2/priscilla/HSparents/scripts/vcf_compare.sh &

#makes a file called all_genotypes_compared.txt with % discordance

#42_20 and 42_23 were not named correctly in subsetted vcf: redo and rerun

vcf-subset --exclude-ref -c 42_23 inbred_68.ug.filter.vsqr.recode.noREP.noINDEL_renamed.vcf > 42_23.inbred.vcf
vcf-subset --exclude-ref -c 42_20 inbred_68.ug.filter.vsqr.recode.noREP.noINDEL_renamed.vcf > 42_20.inbred.vcf


#what is similarity of two totally random populations?
#need to first rename the vcf samples so they have matching names to re-run

bcftools reheader -s test.txt 12BME10_229.inbred.vcf > test.inbred.vcf
bcftools reheader -s test.txt 13_34.hsparent.vcf > test.hsparent.vcf

/mnt/internal_2/priscilla/HSparents/Programs/vcflib/bin/vcfannotategenotypes inbred test.hsparent.vcf test.inbred.vcf > test.annotated.vcf

/mnt/internal_2/priscilla/HSparents/Programs/vcflib/bin/vcfgenotypecompare inbred test.annotated.vcf > test.gtcompare.vcf

/mnt/internal_2/priscilla/HSparents/Programs/vcflib/bin/vcf2tsv -n NA -g test.gtcompare.vcf > test.tsv

Rscript /mnt/internal_2/priscilla/HSparents/scripts/compare_genotypes.Rscript test.tsv
#[1] 0.1803926

bcftools reheader -s test.txt  12LN6_16_B78.inbred.vcf > test.inbred.vcf
bcftools reheader -s test.txt DGRP_28253.hsparent.vcf > test.hsparent.vcf
#[1] 0.03320684

bcftools reheader -s test.txt  4_27.inbred.vcf > test.inbred.vcf
bcftools reheader -s test.txt 12LN6_41_B47.hsparent.vcf > test.hsparent.vcf
#[1] 0.3982004

bcftools reheader -s test.txt  I16.inbred.vcf > test.inbred.vcf
bcftools reheader -s test.txt 20_28.hsparent.vcf > test.hsparent.vcf
#[1] 0.3910519

bcftools reheader -s test.txt  12LN10_32.inbred.vcf > test.inbred.vcf
bcftools reheader -s test.txt DGRP_28253.hsparent.vcf > test.hsparent.vcf
[1] 0.07533122

#Wrote a script to conduct simulations that randomly pick files and do the comparison.

nohup /mnt/internal_2/priscilla/HSparents/scripts/vcf_compare_simulation.sh &

# update 7/7/17. I just realized that I had an error in my barcode-sample sheet. Wells F1-F4 were mislabeled. The samples affected are 4_12 (which was not used in the swarm), 4_27, 13_29, and 24_9. Also, 24_2 was omitted.

#first will remove the incorrect bam files from the mapped folder
cd /mnt/internal_2/priscilla/HSparents/mapped
rm *24_9*
rm *4_12*
rm *13_29*
rm *4_27*

#also remove the 4_12 files from the inbred_mapped folder
cd ../inbred_mapped
rm *4_12*

#move all files in the variants folder to a subdirectory called “old”
cd ../variants
mkdir old
mv * old

#now can rework all variants files without any old files interfering

#make a file in RawData called samples_redo.txt that has the fastq filenames that were affected
index_8nt_682_SL254590.fastq.gz
index_8nt_685_SL254591.fastq.gz
1978_ARL_687_SL254592.fastq.gz
1978_ARL_691_SL254593.fastq.gz

#remake a corrected index file and sort
awk '{print $1,$2,$3}' hybrid_swarm_index_sample_corrected_aq.txt | sort -k2 > index_sorted_HKCHFALXX.txt

#join together file info and sample info into one file (all_file_info_txt)
join -1 2 -2 2 filenames_HKCHFALXX_sorted.txt index_sorted_HKCHFALXX.txt > HKCHFALXX_all_file_info.txt

#edit the pear_bwa script to map only on the redo samples
#run script that includes parallel command to do all mapping of problems
/mnt/internal_2/priscilla/HSparents/scripts/pear_bwa.sh
ls /mnt/internal_2/priscilla/HSparents/mapped/*.sort.dedup.bam > /mnt/internal_2/priscilla/HSparents/scripts/bams.list

nohup java -jar /usr/local/bin/GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-nt 4 \
-nct 6 \
-hets 0.01 \
-R /mnt/internal_2/priscilla/ref_files/all_dmel.fasta \
-I /mnt/internal_2/priscilla/HSparents/scripts/bams.list  \
-o /mnt/internal_2/priscilla/HSparents/variants/hsparents.ug.raw.vcf &

#[1] 16492

 ### run VariantRecalibrator
nohup java -jar /usr/local/bin/GenomeAnalysisTK.jar \
     -T VariantRecalibrator \
-R /mnt/internal_2/priscilla/ref_files/all_dmel.fasta \
-input /mnt/internal_2/priscilla/HSparents/variants/hsparents.ug.raw.vcf \
-resource:dgrp,known=false,training=true,truth=true,prior=15.0 /mnt/internal_2/priscilla/HSparents/variants/inbred.ug.dgrp2_snp.randomSubset.noREP.noINDEL.vcf.truthiness.vcf \
-an DP \
-an QD \
-an FS \
-an SOR \
-an MQ \
-an MQRankSum \
-an ReadPosRankSum \
-mode SNP \
-nt 10 \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
-recalFile /mnt/internal_2/priscilla/HSparents/variants/recalibrate_SNP.recal \
-tranchesFile /mnt/internal_2/priscilla/HSparents/variants/recalibrate_SNP.tranches \
-rscriptFile /mnt/internal_2/priscilla/HSparents/variants/recalibrate_SNP_plots.R &


### run ApplyRecalibration
nohup java -jar /usr/local/bin/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R /mnt/internal_2/priscilla/ref_files/all_dmel.fasta \
-input /mnt/internal_2/priscilla/HSparents/variants/hsparents.ug.raw.vcf \
-mode SNP \
--ts_filter_level 99.0 \
-recalFile /mnt/internal_2/priscilla/HSparents/variants/recalibrate_SNP.recal \
-tranchesFile /mnt/internal_2/priscilla/HSparents/variants/recalibrate_SNP.tranches \
-o /mnt/internal_2/priscilla/HSparents/variants/hsparents.ug.raw.vsqr.vcf &

### filter sites based on PASS field
vcftools \
--vcf /mnt/internal_2/priscilla/HSparents/variants/hsparents.ug.raw.vsqr.vcf \
--remove-filtered-all \
--recode \
--recode-INFO-all \
--out /mnt/internal_2/priscilla/HSparents/variants/hsparents.ug.filter.vsqr

### remove repetitive regions & filter around indels & convert to BCF
bedtools intersect -sorted -v -header \
-a /mnt/internal_2/priscilla/HSparents/variants/hsparents.ug.filter.vsqr.recode.vcf \
-b /mnt/internal_2/priscilla/HSparents/variants/repMasker.dgrp2_indel_100bp.sort.bed | \
tee /mnt/internal_2/priscilla/HSparents/variants/hsparents.ug.filter.vsqr.recode.noREP.noINDEL.vcf | \
bgzip -c > \
/mnt/internal_2/priscilla/HSparents/variants/hsparents.filter.vsqr.recode.noREP.noINDEL.vcf.gz

ls *merge.bam > ../scripts/bam1.txt
ls *up.bam > ../scripts/bam2.txt
cat bam1.txt bam2.txt > inbred_bams.list



nohup java -jar /usr/local/bin/GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-nt 4 \
-nct 6 \
-hets 0.01 \
-R /mnt/internal_2/priscilla/ref_files/all_dmel.fasta \
-I /mnt/internal_2/priscilla/HSparents/scripts/inbred_bams.list  \
-o /mnt/internal_2/priscilla/HSparents/variants/inbred_68.ug.raw.vcf &


#combine remaining processes into vcf_inbred.sh to run as one
 ### run VariantRecalibrator
nohup java -jar /usr/local/bin/GenomeAnalysisTK.jar \
     -T VariantRecalibrator \
-R /mnt/internal_2/priscilla/ref_files/all_dmel.fasta \
-input /mnt/internal_2/priscilla/HSparents/variants/inbred_68.ug.raw.vcf \
-resource:dgrp,known=false,training=true,truth=true,prior=15.0 /mnt/internal_2/priscilla/HSparents/variants/inbred.ug.dgrp2_snp.randomSubset.noREP.noINDEL.vcf.truthiness.vcf \
-an DP \
-an QD \
-an FS \
-an SOR \
-an MQ \
-an MQRankSum \
-an ReadPosRankSum \
-mode SNP \
-nt 10 \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
-recalFile /mnt/internal_2/priscilla/HSparents/variants/recalibrate_SNP_inbred_68.recal \
-tranchesFile /mnt/internal_2/priscilla/HSparents/variants/recalibrate_SNP.tranches \
-rscriptFile /mnt/internal_2/priscilla/HSparents/variants/recalibrate_SNP_plots_inbred_68.R &

nohup java -jar /usr/local/bin/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R /mnt/internal_2/priscilla/ref_files/all_dmel.fasta \
-input /mnt/internal_2/priscilla/HSparents/variants/inbred_68.ug.raw.vcf \
-mode SNP \
--ts_filter_level 99.0 \
-recalFile /mnt/internal_2/priscilla/HSparents/variants/recalibrate_SNP_inbred_68.recal \
-tranchesFile /mnt/internal_2/priscilla/HSparents/variants/recalibrate_SNP.tranches \
-o /mnt/internal_2/priscilla/HSparents/variants/inbred_68.ug.raw.vsqr.vcf &

### filter sites based on PASS field
vcftools \
--vcf /mnt/internal_2/priscilla/HSparents/variants/inbred_68.ug.raw.vsqr.vcf \
--remove-filtered-all \
--recode \
--recode-INFO-all \
--out /mnt/internal_2/priscilla/HSparents/variants/inbred_68.ug.filter.vsqr

### remove repetitive regions & filter around indels & convert to BCF
bedtools intersect -sorted -v -header \
-a /mnt/internal_2/priscilla/HSparents/variants/inbred_68.ug.filter.vsqr.recode.vcf \
-b /mnt/internal_2/priscilla/HSparents/variants/repMasker.dgrp2_indel_100bp.sort.bed | \
tee /mnt/internal_2/priscilla/HSparents/variants/inbred_68.ug.filter.vsqr.recode.noREP.noINDEL.vcf | \
bgzip -c > \
/mnt/internal_2/priscilla/HSparents/variants/inbred_68.ug.filter.vsqr.recode.noREP.noINDEL.vcf.gz

#rename inbred vcf sample names and process both files into individual files 
bcftools query -l inbred_68.ug.filter.vsqr.recode.noREP.noINDEL.vcf.gz >inbred_68_samples.txt

cp inbred_68_samples.txt inbred_68_samples_edit.txt

sed -i 's/Ithica_//g' inbred_68_samples_edit.txt
sed -i 's/\-/\_/g' inbred_68_samples_edit.txt
sed -i 's/(//g' inbred_68_samples_edit.txt
sed -i 's/)//g' inbred_68_samples_edit.txt
sed -i 's/Birmingham_AL_//g' inbred_68_samples_edit.txt
sed -i 's/Selba_AL_//g' inbred_68_samples_edit.txt
sed -i 's/TampaBay_FL_//g' inbred_68_samples_edit.txt
sed -i 's/Thomasville_GA_//g' inbred_68_samples_edit.txt
sed -i 's/BullocksHarbor_BerryIslands_//g' inbred_68_samples_edit.txt
sed -i 's/Cockburn_SanSalvador_//g' inbred_68_samples_edit.txt
sed -i 's/Freeport_GrandBahamasWest_//g' inbred_68_samples_edit.txt
sed -i 's/GeorgeTown_Exumas_//g' inbred_68_samples_edit.txt
sed -i 's/Mayaguana_Mayaguana_//g' inbred_68_samples_edit.txt
sed -i 's/Meridian_MS_//g' inbred_68_samples_edit.txt

#make copy of inbred_68 vcf file
cp inbred_68.ug.filter.vsqr.recode.noREP.noINDEL.vcf.gz inbred_68.ug.filter.vsqr.recode.noREP.noINDEL_copy.vcf.gz

#attempt to rename samples in bcftools

bcftools reheader -s inbred_68_samples_edit.txt inbred_68.ug.filter.vsqr.recode.noREP.noINDEL_copy.vcf.gz > inbred_68.ug.filter.vsqr.recode.noREP.noINDEL_renamed.vcf.gz

#unzip and keep a zipped copy
gunzip -k inbred_68.ug.filter.vsqr.recode.noREP.noINDEL_renamed.vcf.gz

#Split each vcf into separate vcfs to use gtcompare tool;


bcftools query -l inbred_68.ug.filter.vsqr.recode.noREP.noINDEL_renamed.vcf | parallel -j 8 'vcf-subset --exclude-ref -c {} inbred_68.ug.filter.vsqr.recode.noREP.noINDEL_renamed.vcf > {}.inbred.vcf'


bcftools query -l hsparents.ug.filter.vsqr.recode.noREP.noINDEL.vcf | parallel -j 8 'vcf-subset --exclude-ref -c {} hsparents.ug.filter.vsqr.recode.noREP.noINDEL.vcf > {}.hsparent.vcf'

#mistake in names of 42_@0 and 42_23. Edit sample names file, rename vcf and redo split of these two
vcf-subset --exclude-ref -c 42_20 inbred_68.ug.filter.vsqr.recode.noREP.noINDEL_renamed.vcf > 42_20.inbred.vcf

vcf-subset --exclude-ref -c 42_23 inbred_68.ug.filter.vsqr.recode.noREP.noINDEL_renamed.vcf > 42_23.inbred.vcf

#NEED TO NOT USE --exclude-ref in vcf-split command!!

#Split each vcf into separate vcfs to use gtcompare tool;


bcftools query -l inbred_68.ug.filter.vsqr.recode.noREP.noINDEL_renamed.vcf | parallel -j 8 'vcf-subset -c {} inbred_68.ug.filter.vsqr.recode.noREP.noINDEL_renamed.vcf > {}.inbred.vcf'


bcftools query -l hsparents.ug.filter.vsqr.recode.noREP.noINDEL.vcf | parallel -j 8 'vcf-subset -c {} hsparents.ug.filter.vsqr.recode.noREP.noINDEL.vcf > {}.hsparent.vcf'


#find out if 20_17 and 20_28 are switched using program: vcf_compare_different.sh

#no, they are not switched. They are more similar to the expected files than to the opposite, so maybe one data set is just not great quality for them.

#get list of 12LN6 samples in Alan’s vcf:
bcftools query -l inbred.ug.filter.vsqr.recode.noREP.noINDEL.vcf | grep 12LN6 >12LN6_files_to_test.txt

bcftools query -l inbred.ug.filter.vsqr.recode.noREP.noINDEL.vcf | grep 12LN6 | parallel -j 8 'vcf-subset -c {} inbred.ug.filter.vsqr.recode.noREP.noINDEL.vcf > {}.inbred.vcf'

#one line of Alan’s vcf has a huge string of DNA sequence in it that is messing things up. can’t seem to work around it, so maybe remake a vcf of just the 12LN6 lines to do the comparison?


ls /mnt/icy_1/inbredLines/01_mapped/12LN6.12LN6*.rmdup.bam > /mnt/internal_2/priscilla/HSparents/scripts/12LN6_bams.list

cd /mnt/icy_1/inbredLines/01_mapped
nohup java -jar /usr/local/bin/GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-nt 4 \
-nct 6 \
-hets 0.01 \
-R /mnt/internal_2/priscilla/ref_files/all_dmel.fasta \
-I /mnt/internal_2/priscilla/HSparents/scripts/12LN6_bams.list  \
-o /mnt/internal_2/priscilla/HSparents/variants/12LN6.ug.raw.vcf &

#[1] 11912

#make vcf_12LN6.sh file to finish:

 ### run VariantRecalibrator
nohup java -jar /usr/local/bin/GenomeAnalysisTK.jar \
     -T VariantRecalibrator \
-R /mnt/internal_2/priscilla/ref_files/all_dmel.fasta \
-input /mnt/internal_2/priscilla/HSparents/variants/12LN6.ug.raw.vcf \
-resource:dgrp,known=false,training=true,truth=true,prior=15.0 /mnt/internal_2/priscilla/HSparents/variants/inbred.ug.dgrp2_snp.randomSubset.noREP.noINDEL.vcf.truthiness.vcf \
-an DP \
-an QD \
-an FS \
-an SOR \
-an MQ \
-an MQRankSum \
-an ReadPosRankSum \
-mode SNP \
-nt 10 \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
-recalFile /mnt/internal_2/priscilla/HSparents/variants/recalibrate_SNP_12LN6.recal \
-tranchesFile /mnt/internal_2/priscilla/HSparents/variants/recalibrate_SNP.tranches \
-rscriptFile /mnt/internal_2/priscilla/HSparents/variants/recalibrate_SNP_plots_12LN6.R &

nohup java -jar /usr/local/bin/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R /mnt/internal_2/priscilla/ref_files/all_dmel.fasta \
-input /mnt/internal_2/priscilla/HSparents/variants/12LN6.ug.raw.vcf \
-mode SNP \
--ts_filter_level 99.0 \
-recalFile /mnt/internal_2/priscilla/HSparents/variants/recalibrate_SNP_12LN6.recal \
-tranchesFile /mnt/internal_2/priscilla/HSparents/variants/recalibrate_SNP.tranches \
-o /mnt/internal_2/priscilla/HSparents/variants/12LN6.ug.raw.vsqr.vcf &

### filter sites based on PASS field
vcftools \
--vcf /mnt/internal_2/priscilla/HSparents/variants/12LN6.ug.raw.vsqr.vcf \
--remove-filtered-all \
--recode \
--recode-INFO-all \
--out /mnt/internal_2/priscilla/HSparents/variants/12LN6.ug.filter.vsqr

### remove repetitive regions & filter around indels & convert to BCF
bedtools intersect -sorted -v -header \
-a /mnt/internal_2/priscilla/HSparents/variants/12LN6.ug.filter.vsqr.recode.vcf \
-b /mnt/internal_2/priscilla/HSparents/variants/repMasker.dgrp2_indel_100bp.sort.bed | \
tee /mnt/internal_2/priscilla/HSparents/variants/12LN6.ug.filter.vsqr.recode.noREP.noINDEL.vcf | \
bgzip -c > \
/mnt/internal_2/priscilla/HSparents/variants/12LN6.ug.filter.vsqr.recode.noREP.noINDEL.vcf.gz

bcftools query -l 12LN6.ug.filter.vsqr.recode.noREP.noINDEL.vcf > 12LN6_sample_names.txt

sed -i 's/\-/\_/g' 12LN6_sample_names.txt
sed -i 's/(//g' 12LN6_sample_names.txt
sed -i 's/)//g' 12LN6_sample_names.txt

bcftools reheader -s 12LN6_sample_names.txt 12LN6.ug.filter.vsqr.recode.noREP.noINDEL.vcf.gz > 12LN6.ug.filter.vsqr.recode.noREP.noINDEL.renamed.vcf.gz

bcftools query -l 12LN6.ug.filter.vsqr.recode.noREP.noINDEL.renamed.vcf.gz | parallel -j 8 'vcf-subset -c {} 12LN6.ug.filter.vsqr.recode.noREP.noINDEL.renamed.vcf.gz > /mnt/internal_2/priscilla/HSparents/variants/12LN6/{}.inbred.vcf'



###########

#9/13/17 
#Need to change read group info in all of the “known” inbred files so that they can be merged with mine.

#List of all files can be found in /mnt/internal_2/priscilla/HSparents/scripts/inbred_bams.list
#List of old names and new names can be found at 
