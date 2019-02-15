#!/bin/bash

masterDir="/scratch/pae3g/hybrid"
echo masterDir = $masterDir

GATK="$masterDir/etc/GenomeAnalysisTK.jar"
echo GATK = $GATK

inputFileStem=${1}
inputDir="/scratch/pae3g/parents/bams"

output_gVCF="/scratch/pae3g/parents/gVCF"

referenceGenome="$masterDir""/seqs/all_dmel.fasta"

    java -jar -Xmx2G $GATK \
		     -T HaplotypeCaller \
		     -R $referenceGenome \
		     -I $inputDir/${inputFileStem}.bams.list \
		     -hets 0.01 \
		     -indelHeterozygosity 0.001 \
		     --emitRefConfidence GVCF \
		     --variant_index_type LINEAR \
		     --variant_index_parameter 128000 \
		     -o $output_gVCF/${inputFileStem}.raw.snps.indels.g.vcf



