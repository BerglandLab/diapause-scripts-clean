#!/bin/bash

# let us stop the whole script on CTRL+C
	int_handler()
	{
		echo "Interrupted."
		# Kill the parent process of the script.
		kill $PPID
		exit 1
	}
	trap 'int_handler' INT

### functions

	
	mapReads () {

		## define some parameters
			sample=${1}			
			run=${2}
			inputDir=${3}
		
			threads=10
		
			outputDir="/mnt/inbred/inbredLines/01_mapped"
	
			lanes=($( ls ${inputDir}/${run}/FASTQ/*${sample}* | grep -oE 'L[0-9]{1,}' | sort | uniq ))
	
			strain=$(grep ${sample} ${inputDir}/${run}/*.csv | cut -d',' -f1)
			
				
		## map reads
			for lane in "${lanes[@]}"; do
		
				### first, merge overlapping reads
					pear-0.9.6-bin-64 \
					-f ${inputDir}/${run}/FASTQ/*${sample}*${lane}_R1_* \
					-r ${inputDir}/${run}/FASTQ/*${sample}*${lane}_R2_* \
					-o ${inputDir}/${run}/FASTQ/${sample}_${lane} \
					-j ${threads}
		
				### next, map to reference genome
					## assembled reads (i.e., those that were merged by PEAR)
						bwa mem -t ${threads} \
							-R "@RG\tID:${strain};${run};${lane};unassembled\tSM:${strain}" \
							/mnt/Alan/dmel_reference/all_dmel.fasta \
							 ${inputDir}/${run}/FASTQ/${sample}_${lane}.assembled.fastq | \
						samtools view -Suh - | \
						samtools sort -@ ${threads} - ${outputDir}/${strain}.${run}.${lane}.assembled.sort
				
						samtools index ${outputDir}/${strain}.${run}.${lane}.assembled.sort.bam
						/mnt/inbred/inbredLines/software/samtools-0.1.19/samtools rmdup -s ${outputDir}/${strain}.${run}.${lane}.assembled.sort.bam ${outputDir}/${strain}.${run}.${lane}.assembled.sort.rmdup.bam
						samtools index ${outputDir}/${strain}.${run}.${lane}.assembled.sort.rmdup.bam
						
					## unassembled reads
						bwa mem -t ${threads} \
							-R "@RG\tID:${strain};${run};${lane};assembled\tSM:${strain}" \
							/mnt/Alan/dmel_reference/all_dmel.fasta \
							 ${inputDir}/${run}/FASTQ/${sample}_${lane}.unassembled.forward.fastq \
							 ${inputDir}/${run}/FASTQ/${sample}_${lane}.unassembled.reverse.fastq | \
						samtools view -Suh - | \
						samtools sort -@ ${threads} - ${outputDir}/${strain}.${run}.${lane}.unassembled.sort
				
						samtools index ${outputDir}/${strain}.${run}.${lane}.unassembled.sort.bam
						/mnt/inbred/inbredLines/software/samtools-0.1.19/samtools rmdup ${outputDir}/${strain}.${run}.${lane}.unassembled.sort.bam ${outputDir}/${strain}.${run}.${lane}.unassembled.sort.rmdup.bam
						samtools index ${outputDir}/${strain}.${run}.${lane}.unassembled.sort.rmdup.bam
					
						
			done
	
	
	}
	export -f mapReads



	mapReads_sraData () {

		## define some parameters
			dir=${1}
			pop=($( echo $dir | cut -f 6 -d "/" ))		
			strains=($( ls $dir ))
		
			threads=10
		
			outputDir="/mnt/inbred/inbredLines/01_mapped"
	
				
		## map reads
			for strain in "${strains[@]}"; do

				runs=($(ls ${dir}/${strain}/*.gz | grep -Eo 'SRR[0-9]{1,}' | sort | uniq))
				
				for run in "${runs[@]}"; do

					files=($(ls ${dir}/${strain}/${run}*.gz))
					
					if [ ${#files[@]} = "1" ]; then
					
						format=$(gzip -dc ${dir}/${strain}/${run}_1.fastq.gz | head | /mnt/inbred/inbredLines/software/DetermineFastqQualityEncoding.pl - 2> foo | grep "format" | grep -Eo "64|33")
						
						if [ ${format} = "64" ]; then
							gzip -dc ${dir}/${strain}/${run}_1.fastq.gz | /mnt/inbred/inbredLines/software/convFQ.perl | gzip > ${dir}/${strain}/${run}_1.conv.fastq.gz
							rm 	${dir}/${strain}/${run}_1.fastq.gz
							mv ${dir}/${strain}/${run}_1.conv.fastq.gz ${dir}/${strain}/${run}_1.fastq.gz
						fi
							
						bwa mem -t ${threads} \
						-R "@RG\tID:${pop}_${strain};${run};se\tSM:${pop}_${strain}" \
						/mnt/Alan/dmel_reference/all_dmel.fasta \
						${dir}/${strain}/${run}_1.fastq.gz | \
						samtools view -Suh - | \
						samtools sort -@ ${threads} - ${outputDir}/${pop}.${strain}.${run}.sort
						
			
						samtools index ${outputDir}/${pop}.${strain}.${run}.sort.bam
						/mnt/inbred/inbredLines/software/samtools-0.1.19/samtools rmdup -s ${outputDir}/${pop}.${strain}.${run}.sort.bam ${outputDir}/${pop}.${strain}.${run}.sort.rmdup.bam
						samtools index ${outputDir}/${pop}.${strain}.${run}.sort.rmdup.bam
						

					else
					
						format=$(gzip -dc ${dir}/${strain}/${run}_1.fastq.gz | head | /mnt/inbred/inbredLines/software/DetermineFastqQualityEncoding.pl - 2> foo | grep "format" | grep -Eo "64|33")
						if [ ${format} = "64" ]; then
							gzip -dc ${dir}/${strain}/${run}_1.fastq.gz | /mnt/inbred/inbredLines/software/convFQ.perl | gzip > ${dir}/${strain}/${run}_1.conv.fastq.gz
							rm 	${dir}/${strain}/${run}_1.fastq.gz
							mv ${dir}/${strain}/${run}_1.conv.fastq.gz ${dir}/${strain}/${run}_1.fastq.gz
							
							gzip -dc ${dir}/${strain}/${run}_2.fastq.gz | /mnt/inbred/inbredLines/software/convFQ.perl | gzip > ${dir}/${strain}/${run}_2.conv.fastq.gz
							rm 	${dir}/${strain}/${run}_2.fastq.gz
							mv ${dir}/${strain}/${run}_2.conv.fastq.gz ${dir}/${strain}/${run}_2.fastq.gz

						fi
					
					
						bwa mem -t ${threads} \
						-R "@RG\tID:${pop}_${strain};${run};pe\tSM:${pop}_${strain}" \
						/mnt/Alan/dmel_reference/all_dmel.fasta \
						${dir}/${strain}/${run}_1.fastq.gz \
						${dir}/${strain}/${run}_2.fastq.gz | \
						samtools view -Suh - | \
						samtools sort -@ ${threads} - ${outputDir}/${pop}.${strain}.${run}.sort
						
						samtools index ${outputDir}/${pop}.${strain}.${run}.sort.bam
						/mnt/inbred/inbredLines/software/samtools-0.1.19/samtools rmdup ${outputDir}/${pop}.${strain}.${run}.sort.bam ${outputDir}/${pop}.${strain}.${run}.sort.rmdup.bam
						samtools index ${outputDir}/${pop}.${strain}.${run}.sort.rmdup.bam

					fi
					 
				done
							
			done
	
	
	}
	export -f mapReads_sraData





	parseVCF () {
		
		file=${1}
		
		cat ${file} | awk -v file=${file} '{
			if($1=="#CHROM") {
				printf $1","$2",ref,alt," > file".allele.delim"
				printf $1","$2",ref,alt," > file".rd.delim"
				
				for(i=10; i<=NF; i++) {
					printf $i > file".allele.delim"
					printf $i > file".rd.delim"
					
					if(i<NF) printf "," > file".allele.delim"
					if(i<NF) printf "," > file".rd.delim"
					
					if(i==NF) printf "\n" > file".allele.delim"
					if(i==NF) printf "\n" > file".rd.delim"

				}
			} else if (substr($1, 1, 1)!="#" && length($5)==1) {
				printf $1","$2","$4","$5","  > file".allele.delim"
				printf $1","$2","$4","$5","  > file".rd.delim"
				
				for(i=10; i<=NF; i++) {
					split($i, sp, ":")

					## alleles 
						if(sp[1]=="0/0") printf "0" > file".allele.delim"
						if(sp[1]=="1/1") printf "1" > file".allele.delim"
						if(sp[1]=="1/0" || sp[1]=="0/1") printf "0.5" > file".allele.delim"
						if(sp[1]=="./.") printf "NA" > file".allele.delim"
					
						if(i<NF) printf "," > file".allele.delim"
						if(i==NF) printf "\n" > file".allele.delim"
						
					## read depth
						if(sp[1]!="./.") {
							printf sp[3] > file".rd.delim"
						} else {
							printf 0 > file".rd.delim"
						}
						if(i<NF) printf "," > file".rd.delim"
						if(i==NF) printf "\n" > file".rd.delim"
						
				}
					
			}
		}'
	}
	export -f parseVCF
	
	renameBam () {
		origFile=$(echo ${1} | cut -d'/' -f 6)
		pop=${2}
		mv /mnt/inbred/inbredLines/01_mapped/${origFile} /mnt/inbred/inbredLines/01_mapped/${pop}.${origFile}
	}
	export -f renameBam

			
	gatk_HC_pipeline () {

		inputDir=/mnt/inbred/inbredLines/01_mapped
		inputFileStem=${1}
		output_gVCF=/mnt/inbred/inbredLines/02_gVCF
		
		ls ${inputDir}/${inputFileStem}*.rmdup.bam > /mnt/inbred/inbredLines/scripts/${inputFileStem}_bams.list
		
		cat /mnt/inbred/inbredLines/scripts/${inputFileStem}_bams.list | awk '{
			n=split($1, stem, "\\.")
			split($1, file, "\\/")
			
			printf file[6]"\t"
			for(i=1; i<n; i++) {
				printf stem[i]"."
			} 
			printf "realigned.bam"
			printf "\n"
		}' > /mnt/inbred/inbredLines/scripts/${inputFileStem}_input_output.map
			
		java -Xmx1G -jar ~/GenomeAnalysisTK.jar  \
		-T RealignerTargetCreator \
		-R /mnt/Alan/dmel_reference/all_dmel.fasta \
		-I /mnt/inbred/inbredLines/scripts/${inputFileStem}_bams.list \
		-o ${inputDir}/${inputFileStem}.intervals
			
		java -Xmx1G -jar ~/GenomeAnalysisTK.jar \
		-T IndelRealigner \
		-R /mnt/Alan/dmel_reference/all_dmel.fasta \
		-I /mnt/inbred/inbredLines/scripts/${inputFileStem}_bams.list \
		--targetIntervals ${inputDir}/${inputFileStem}.intervals \
		--nWayOut /mnt/inbred/inbredLines/scripts/${inputFileStem}_input_output.map
						
		# rm ${inputDir}/${inputFileStem}*.rmdup.bam
		ls ${inputDir}/${inputFileStem}*.rmdup.realigned.bam > /mnt/inbred/inbredLines/scripts/${inputFileStem}_bams.list
			
		java -Xmx1G -jar ~/GenomeAnalysisTK.jar \
		-T HaplotypeCaller \
		-R /mnt/Alan/dmel_reference/all_dmel.fasta \
		-I /mnt/inbred/inbredLines/scripts/${inputFileStem}_bams.list \
		-hets 0.01 \
		-indelHeterozygosity 0.001 \
		--emitRefConfidence GVCF \
		--variant_index_type LINEAR \
		--variant_index_parameter 128000 \
		-o /mnt/inbred/inbredLines/02_gVCF/${inputFileStem}.rmdup.realign.raw.snps.indels.g.vcf
   
   }
   export -f gatk_HC_pipeline     
        
   runfreebayes () {
	
		ulimit -n 10000
	
		chr=$(echo ${1} | cut -f 1 -d":")
	
		freebayes \
		-L /mnt/inbred/inbredLines/scripts/bams.list \
		-v /mnt/inbred/inbredLines/02_variants_freebayes/inbred.${chr}.vcf \
		-f /mnt/Alan/dmel_reference/all_dmel.fasta \
		--populations /mnt/inbred/inbredLines/scripts/pops.list \
		-T 0.01 \
		-r ${1} \
		-Q 20 \
		-m 10 \
		-w
	}
	export -f runfreebayes
		
### pipeline
	
	### first, map all the reads
		### first plate, first run		
			# samps=($(cat /mnt/inbred/inbredLines/00_rawData/ab_03_17_2015/ab_03_17_2015.plate1.info.csv | cut -d',' -f8 | grep -v "NGX"))
			# parallel --gnu -j1 mapReads ::: "${samps[@]}" ::: "ab_03_17_2015" ::: "/mnt/inbred/inbredLines/00_rawData"
		
		### second plate, first run
			#samps=($(cat /mnt/inbred/inbredLines/00_rawData/plate_2/plate2_inbreedingstrains.info.csv | cut -d',' -f8 | grep -v "NGX"))
			#parallel --gnu -j1 mapReads ::: "${samps[@]}" ::: "plate_2" ::: "/mnt/inbred/inbredLines/00_rawData"
		
		### third plate, first run
			# samps=($(cat /mnt/inbred/inbredLines/00_rawData/plate_3/Plate3_Inbreeding_Strains.info.csv | cut -d',' -f8 | grep -v "NGX"))
			# parallel --gnu -j1 mapReads ::: "${samps[@]}" ::: "plate_3" ::: "/mnt/inbred/inbredLines/00_rawData"

		### fourth plate, first run
			# samps=($(cat /mnt/inbred/inbredLines/00_rawData/plate_4/Plate_4_Inbreeding_Strains.csv | cut -d',' -f8 | grep -v "NGX"))
			# parallel --gnu -j1 mapReads ::: "${samps[@]}" ::: "plate_4" ::: "/mnt/inbred/inbredLines/00_rawData"

		### rename ME & PA bam files
			# parallel --gnu -j1 renameBam ::: $(ls /mnt/inbred/inbredLines/01_mapped/12BME10*) ::: 12BME10
			# parallel --gnu -j1 renameBam ::: $(ls /mnt/inbred/inbredLines/01_mapped/12LN6*) ::: 12LN6
			# parallel --gnu -j1 renameBam ::: $(ls /mnt/inbred/inbredLines/01_mapped/LNPA*) ::: 12LN6
			# parallel --gnu -j1 renameBam ::: $(ls /mnt/inbred/inbredLines/01_mapped/12LN10*) ::: 12LN10
			# parallel --gnu -j1 renameBam ::: $(ls /mnt/inbred/inbredLines/01_mapped/12LNPA*) ::: 12LN6
			# parallel --gnu -j1 renameBam ::: $(ls /mnt/inbred/inbredLines/01_mapped/\(Unknown\)*) ::: UnkNAM

		### external data
			# parallel --gnu -j1 mapReads_sraData ::: "/mnt/inbred/inbredLines/00_rawData/Birmingham_AL"
			# parallel --gnu -j1 mapReads_sraData ::: "/mnt/inbred/inbredLines/00_rawData/BullocksHarbor_BerryIslands"
			# parallel --gnu -j1 mapReads_sraData ::: "/mnt/inbred/inbredLines/00_rawData/Cockburn_SanSalvador"
			# parallel --gnu -j1 mapReads_sraData ::: "/mnt/inbred/inbredLines/00_rawData/dgrp"
			# parallel --gnu -j1 mapReads_sraData ::: "/mnt/inbred/inbredLines/00_rawData/Freeport_GrandBahamasWest"
			# parallel --gnu -j1 mapReads_sraData ::: "/mnt/inbred/inbredLines/00_rawData/GeorgeTown_Exumas"
			# parallel --gnu -j1 mapReads_sraData ::: "/mnt/inbred/inbredLines/00_rawData/Ithica"
			# parallel --gnu -j1 mapReads_sraData ::: "/mnt/inbred/inbredLines/00_rawData/Mayaguana_Mayaguana"
			# parallel --gnu -j1 mapReads_sraData ::: "/mnt/inbred/inbredLines/00_rawData/Meridian_MS"
			# parallel --gnu -j1 mapReads_sraData ::: "/mnt/inbred/inbredLines/00_rawData/PortAuPrince_Haiti"
			# parallel --gnu -j1 mapReads_sraData ::: "/mnt/inbred/inbredLines/00_rawData/Sebastian_FL"
			# parallel --gnu -j1 mapReads_sraData ::: "/mnt/inbred/inbredLines/00_rawData/Selba_AL"
			# parallel --gnu -j1 mapReads_sraData ::: "/mnt/inbred/inbredLines/00_rawData/TampaBay_FL"
			# parallel --gnu -j1 mapReads_sraData ::: "/mnt/inbred/inbredLines/00_rawData/Thomasville_GA"
			# parallel --gnu -j1 mapReads_sraData ::: "/mnt/inbred/inbredLines/00_rawData/dpgp3"
			# parallel --gnu -j1 mapReads_sraData ::: "/mnt/inbred/inbredLines/00_rawData/CO"
			# parallel --gnu -j1 mapReads_sraData ::: "/mnt/inbred/inbredLines/00_rawData/UG"
			# parallel --gnu -j1 mapReads_sraData ::: "/mnt/inbred/inbredLines/00_rawData/NL"
			# parallel --gnu -j1 mapReads_sraData ::: "/mnt/inbred/inbredLines/00_rawData/France"
			# parallel --gnu -j1 mapReads_sraData ::: "/mnt/inbred/inbredLines/00_rawData/GU"
	
		### next, call variants w/ GATK
#			ls /mnt/inbred/inbredLines/01_mapped/*.sort.rmdup.bam > /mnt/inbred/inbredLines/scripts/bams.list
#	
#			ulimit -n 10000
#				
#			 java -jar ~/GenomeAnalysisTK.jar \
#			 -T UnifiedGenotyper \
#			 -nt 4 \
#			 -nct 6 \
#			 -hets 0.01 \
#			 -R /mnt/Alan/dmel_reference/all_dmel.fasta \
#			 -I /mnt/inbred/inbredLines/scripts/bams.list \
#			 -o /mnt/inbred/inbredLines/02_variants/inbred.ug.raw.vcf
#
#		### run VSQR analysis
#			### make thruthiness set
#				bedtools intersect -sorted -header \
#				-a /mnt/inbred/inbredLines/02_variants/inbred.ug.raw.vcf \
#				-b /mnt/inbred/inbredLines/filterData/dgrp2_snp.randomSubset.noREP.noINDEL.vcf > \
#				/mnt/inbred/inbredLines/02_variants/inbred.ug.dgrp2_snp.randomSubset.noREP.noINDEL.vcf.truthiness.vcf
#
#			### run VariantRecalibrator
#				java -jar ~/GenomeAnalysisTK.jar \
#				-T VariantRecalibrator \
#				-R /mnt/Alan/dmel_reference/all_dmel.fasta \
#				-input /mnt/inbred/inbredLines/02_variants/inbred.ug.raw.vcf \
#				-resource:dgrp,known=false,training=true,truth=true,prior=15.0 /mnt/inbred/inbredLines/02_variants/inbred.ug.dgrp2_snp.randomSubset.noREP.noINDEL.vcf.truthiness.vcf \
#				-an DP \
#				-an QD \
#				-an FS \
#				-an SOR \
#				-an MQ \
#				-an MQRankSum \
#				-an ReadPosRankSum \
#				-mode SNP \
#				-nt 10 \
#				-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
#				-recalFile /mnt/inbred/inbredLines/02_variants/recalibrate_SNP.recal \
#				-tranchesFile /mnt/inbred/inbredLines/02_variants/recalibrate_SNP.tranches \
#				-rscriptFile /mnt/inbred/inbredLines/02_variants/recalibrate_SNP_plots.R
#			
#			### run ApplyRecalibration
#				java -jar ~/GenomeAnalysisTK.jar \
#				-T ApplyRecalibration \
#				-R /mnt/Alan/dmel_reference/all_dmel.fasta \
#				-input /mnt/inbred/inbredLines/02_variants/inbred.ug.raw.vcf \
#				-mode SNP \
#				--ts_filter_level 99.0 \
#				-recalFile /mnt/inbred/inbredLines/02_variants/recalibrate_SNP.recal \
#				-tranchesFile /mnt/inbred/inbredLines/02_variants/recalibrate_SNP.tranches \
#				-o /mnt/inbred/inbredLines/02_variants/inbred.ug.raw.vsqr.vcf
#	
#			### filter sites based on PASS field
#				vcftools \
#				--vcf /mnt/inbred/inbredLines/02_variants/inbred.ug.raw.vsqr.vcf \
#				--remove-filtered-all \
#				--recode \
#				--recode-INFO-all \
#				--out /mnt/inbred/inbredLines/02_variants/inbred.ug.filter.vsqr
#	
			### remove repetative regions & filter around indels & convert to BCF
				bedtools intersect -sorted -v -header \
				-a /mnt/inbred/inbredLines/02_variants/inbred.ug.filter.vsqr.recode.vcf \
				-b /mnt/inbred/inbredLines/filterData/repMasker.dgrp2_indel_100bp.sort.bed | \
				tee /mnt/inbred/inbredLines/02_variants/inbred.ug.filter.vsqr.recode.noREP.noINDEL.vcf | \
				bgzip -c > \
				/mnt/inbred/inbredLines/02_variants/inbred.ug.filter.vsqr.recode.noREP.noINDEL.vcf.gz

			## tabix
				tabix -p vcf /mnt/inbred/inbredLines/02_variants/inbred.ug.filter.vsqr.recode.noREP.noINDEL.vcf.gz
			

	
	
	
	
	
	
	
	
	
	
	
	
	
			
	
	
	
	
	
	
	
	
	
	
	### call variants with FreeBayes
#		ls /mnt/inbred/inbredLines/01_mapped/*.sort.rmdup.bam > /mnt/inbred/inbredLines/scripts/bams.list
#
#		ls /mnt/inbred/inbredLines/01_mapped/*.sort.rmdup.bam | awk '{ 
#			 split($1, sp, "\\/") 
#			 split(sp[6], pop, "\\.") 
#			 print $0"\t" pop[1]
#		}' > /mnt/inbred/inbredLines/scripts/pops.list
#		
#		ulimit -n 10000
#		parallel --gnu -j 9 -a /mnt/inbred/inbredLines/scripts/scafDefn.list runfreebayes 
#	
		

	### next, call variants w/ GATK
		### get list of file stems for each independent sample 
		#	 ls /mnt/inbred/inbredLines/01_mapped/*.sort.rmdup.bam | cut -d'/' -f 6 | cut -d'.' -f1,2 | sort | uniq > /mnt/inbred/inbredLines/scripts/strains.list
	
		### run GATK indel-realignment & HaplotypeCaller for each sample
		#	parallel --gnu -j15 -a /mnt/inbred/inbredLines/scripts/strains.list gatk_HC_pipeline
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	#### this is old Unified Genotyper method
		### next, call variants w/ GATK
#			ulimit -n 10000
#				
#			 java -jar ~/GenomeAnalysisTK.jar \
#			 -T UnifiedGenotyper \
#			 -nt 4 \
#			 -nct 6 \
#			 -R /mnt/Alan/dmel_reference/all_dmel.fasta \
#			 -I /mnt/inbred/inbredLines/02_variants/bams.list \
#			 -o /mnt/inbred/inbredLines/02_variants/inbred.vcf
#
#		### parse vcf
#			parseVCF /mnt/inbred/inbredLines/02_variants/inbred.vcf
#		
#		### make tabixd version of parsed genotype file
#			
#			cat /mnt/inbred/inbredLines/02_variants/inbred.vcf.allele.delim | \
#			awk -F',' '{print $1"\t"$2"\t"$0}'  > /mnt/inbred/inbredLines/02_variants/inbred.vcf.allele.delim.vcf
#
#			rm /mnt/inbred/inbredLines/02_variants/inbred.vcf.allele.delim.vcf.gz
#
#			bgzip /mnt/inbred/inbredLines/02_variants/inbred.vcf.allele.delim.vcf
#			tabix -f -p vcf /mnt/inbred/inbredLines/02_variants/inbred.vcf.allele.delim.vcf.gz
#			
#			
#			
#			
			
			
			
			