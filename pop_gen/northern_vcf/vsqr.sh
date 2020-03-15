#! /bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=64000
#SBATCH --time=48:00:00
#SBATCH --partition=economy
#SBATCH --account=berglandlab

module load java/1.8.0_45

#run VariantRecalibrator
echo "Variant Recalibrator"
java -jar /scratch/ab5dr/GATK/GenomeAnalysisTK.jar \
-T VariantRecalibrator \
-R /scratch/ab5dr/InbredLines/ref/all_dmel.fasta \
-input /scratch/ab5dr/InbredLines/scripts/BME_LN_LNPA.raw.vcf \
-resource:dgrp,known=false,training=true,truth=true,prior=15.0 /scratch/ab5dr/InbredLines/vsqr/inbred.ug.dgrp2_snp.randomSubset.noREP.noINDEL.vcf.truthiness.vcf \
-an DP \
-an QD \
-an FS \
-an SOR \
-an MQ \
-an MQRankSum \
-an ReadPosRankSum \
-mode SNP \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
-recalFile /scratch/ab5dr/InbredLines/vsqr/recalibrate_SNP.recal \
-tranchesFile /scratch/ab5dr/InbredLines/vsqr/recalibrate_SNP.tranches \
-rscriptFile /scratch/ab5dr/InbredLines/vsqr/recalibrate_SNP_plots.R

# run ApplyRecalibration
echo "Apply Recalibration"
java -jar /scratch/ab5dr/GATK/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R /scratch/ab5dr/InbredLines/ref/all_dmel.fasta \
-input /scratch/ab5dr/InbredLines/scripts/BME_LN_LNPA.raw.vcf \
-mode SNP \
--ts_filter_level 99.0 \
-recalFile /scratch/ab5dr/InbredLines/vsqr/recalibrate_SNP.recal \
-tranchesFile /scratch/ab5dr/InbredLines/vsqr/recalibrate_SNP.tranches \
-o /scratch/ab5dr/InbredLines/vsqr/BME_LN_LNPA.raw.vsqr.vcf
