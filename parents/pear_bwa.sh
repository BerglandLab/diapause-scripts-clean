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

#define a function that will run reads through pear and bwa
	mapReads () {

		## define some parameters
	                sample=${1}
			echo "sample=" $sample
			inputDir=${2}
			threads=10
			outputDir="/mnt/internal_2/priscilla/HSparents/mapped"
			run="HKCHFALXX"
			lane="s7"
			strain=$(grep ${sample} ${inputDir}/RawData/HKCHFALXX_all_file_info.txt | cut -d ' ' -f8)
			echo "strain=" $strain

			echo 'starting PEAR'
			pear \
					-f ${inputDir}/RawData/*s7_1_*${sample}* \
					-r ${inputDir}/RawData/*s7_2_*${sample}* \
					-o ${inputDir}/RawData/${sample} \
					-j ${threads}
			
			echo 'mapping assembled reads with bwa and exporting to samtools'
			bwa mem -t ${threads} \
					-R "@RG\tID:${strain};${run};${lane};unassembled\tSM:${strain}" \
					/mnt/internal_2/priscilla/ref_files/all_dmel.fasta \
					${inputDir}/RawData/${sample}.assembled.fastq | \
					samtools view -Suh - | \
					samtools sort -@ ${threads} -o ${outputDir}/${strain}.assembled.sort.bam
			                echo "samtools sorting completed, now indexing"
					samtools index ${outputDir}/${strain}.assembled.sort.bam
					echo "removing PCR duplicates"
					if [ ! -d $outputDir/dedup-report ]; then
						mkdir $outputDir/dedup-report
						fi
	
					java \
					-Xmx20g \
					-Dsnappy.disable=true \
					-jar /usr/local/bin/picard.jar MarkDuplicates \
					REMOVE_DUPLICATES=true \
					I=$outputDir/${strain}.assembled.sort.bam \
					O=$outputDir/${strain}.assembled.sort.dedup.bam \
					M=$outputDir/dedup-report/${strain}.assembled.txt \
					VALIDATION_STRINGENCY=SILENT
					samtools index $outputDir/${strain}.assembled.sort.dedup.bam
 

			echo 'mapping unassembled reads with bwa'			
					## unassembled reads
						bwa mem -t ${threads} \
							-R "@RG\tID:${strain};${run};${lane};assembled\tSM:${strain}" \
							/mnt/internal_2/priscilla/ref_files/all_dmel.fasta \
							 ${inputDir}/RawData/${sample}.unassembled.forward.fastq \
							 ${inputDir}/RawData/${sample}.unassembled.reverse.fastq | \
						samtools view -Suh - | \
						samtools sort -@ ${threads} -o ${outputDir}/${strain}.unassembled.sort.bam
						samtools index ${outputDir}/${strain}.unassembled.sort.bam
				echo 'removing pcr duplicates unassembled reads'
						java \
						-Xmx20g \
						-Dsnappy.disable=true \
						-jar /usr/local/bin/picard.jar MarkDuplicates \
						REMOVE_DUPLICATES=true \
						I=$outputDir/${strain}.unassembled.sort.bam \
						O=$outputDir/${strain}.unassembled.sort.dedup.bam \
						M=$outputDir/dedup-report/${strain}.unassembled.txt \
						VALIDATION_STRINGENCY=SILENT
						samtools index $outputDir/${strain}.unassembled.sort.dedup.bam
 

	}
	
	export -f mapReads
	
	samps=($(cat /mnt/internal_2/priscilla/HSparents/RawData/samples.txt))
	parallel --gnu -j1 mapReads ::: "${samps[@]}" ::: "/mnt/internal_2/priscilla/HSparents"

