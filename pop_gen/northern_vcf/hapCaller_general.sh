#! /bin/bash


echo 'defining variables'
#directory paths currently for Rivanna 
refDir=/scratch/ab5dr/InbredLines/ref
inputDir=/scratch/ab5dr/InbredLines/originalBam
interDir=/scratch/ab5dr/InbredLines/interBam
gvcfOutDir=/scratch/ab5dr/InbredLines/gvcf
knownDir=/scratch/ab5dr/InbredLines/dgrpKnown
GATKDir=/scratch/ab5dr/GATK
PicDir=/scratch/ab5dr/picard 
samDir=/scratch/ab5dr/samtools/samtools-1.2

#indel hotspot identification file generation
echo 'RealignerTargetCreator'
java -jar ${GATKDir}/GenomeAnalysisTK.jar \
   -T RealignerTargetCreator \
   -nct 1 \
   -R ${refDir}/all_dmel.fasta \
   -I ${inputDir}/$1.bam \
   --known ${knownDir}/dgrp2.vcf \
   -o ${interDir}/$1.intervals		

#actual indel realignment 
echo 'IndelRealigner'
java -jar ${GATKDir}/GenomeAnalysisTK.jar \
   -T IndelRealigner \
   -nct 1 \
   -R ${refDir}/all_dmel.fasta \
   -I ${inputDir}/$1.bam \
   -known ${knownDir}/dgrp2.vcf \
   -targetIntervals ${interDir}/$1.intervals \
   -o ${interDir}/$1.realn.bam

#extracting out stuff for RG names
echo 'slice and dice in prep for picard'
SRA=$(echo $1 | cut -d'.' -f3)
IDname=$(echo $1 | cut -d'.' -f1,2)

#replacing read groups 
echo 'replacing read groups'
java -jar ${PicDir}/picard.jar AddOrReplaceReadGroups \
I=${interDir}/$1.realn.bam \
O=${interDir}/$1.realn.rg.bam \
RGID=${IDname} \
RGLB=${SRA} \
RGPL=illumina \
RGPU=unit1 \
RGSM=${IDname} \
SORT_ORDER=coordinate \
CREATE_INDEX=true

#base recalibration
echo 'BaseRecalibrator'
java -jar ${GATKDir}/GenomeAnalysisTK.jar \
   -T BaseRecalibrator \
   -nct 1 \
   -R ${refDir}/all_dmel.fasta \
   -I ${interDir}/$1.realn.rg.bam \
   -knownSites ${knownDir}/dgrp2.vcf \
   -o ${interDir}/$1_recal_data.table 
      
#print reads step - merges bam file with recalibration table     
echo 'PrintReads'
java -jar ${GATKDir}/GenomeAnalysisTK.jar \
   -T PrintReads \
   -nct 1 \
   -R ${refDir}/all_dmel.fasta \
   -I ${interDir}/$1.realn.rg.bam \
   --BQSR ${interDir}/$1_recal_data.table \
   -o ${interDir}/$1.realn.rg.recal.bam 

#haplotype caller 
echo 'HaplotypeCaller'
java -jar ${GATKDir}/GenomeAnalysisTK.jar \
   -R ${refDir}/all_dmel.fasta \
   -T HaplotypeCaller \
   -I ${interDir}/$1.realn.rg.recal.bam \
   --emitRefConfidence GVCF \
   -nct 4 \
   -hets 0.01 \
   -indelHeterozygosity 0.001 \
   -o ${gvcfOutDir}/$1.realn.rg.recal.g.vcf				

echo 'Removing intermediate files'
rm {interDir}/$1*
