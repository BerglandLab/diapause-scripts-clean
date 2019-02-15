#!/bin/bash
mapReads_sraData () {

		## define some parameters
			dir=${1}		
			strains=($( ls $dir ))
		
			threads=10
		
			outputDir="/mnt/internal_2/priscilla/HSparents/inbred_mapped"
	
				
		## map reads
			for strain in "${strains[@]}"; do
			    echo "strain is " $strain
				runs=($(ls ${dir}/${strain}/*.gz | grep -Eo 'SRR[0-9]{1,}' | sort | uniq))
				echo "runs are " $runs
				for run in "${runs[@]}"; do

					files=($(ls ${dir}/${strain}/${run}*.gz))

					
					if [ ${#files[@]} = "1" ]; then
					echo "working on single end data"
						format=$(gzip -dc ${dir}/${strain}/${run}_pass.fastq.gz | head | /mnt/icy_1/inbredLines/software/DetermineFastqQualityEncoding.pl - 2> foo | grep "format" | grep -Eo "64|33")
						echo "format = " $format
						if [ ${format} = "64" ]; then
							gzip -dc ${dir}/${strain}/${run}_pass.fastq.gz | /mnt/icy_1/inbredLines/software/convFQ.perl | gzip > ${dir}/${strain}/${run}_1.conv.fastq.gz
							rm 	${dir}/${strain}/${run}_pass.fastq.gz
							mv ${dir}/${strain}/${run}_1.conv.fastq.gz ${dir}/${strain}/${run}_pass.fastq.gz
						fi

						echo "mapping single end reads"
						bwa mem -t ${threads} \
						-R "@RG\tID:${strain};${run};se\tSM:${strain}" \
						/mnt/internal_2/priscilla/ref_files/all_dmel.fasta \
						${dir}/${strain}/${run}_pass.fastq.gz | \
						samtools view -Suh - | \
						samtools sort -@ ${threads} -o ${outputDir}/${strain}.${run}.sort.bam
						
						echo "indexing bam"
						samtools index ${outputDir}/${strain}.${run}.sort.bam
					
						echo "removing PCR duplicates"
						if [ ! -d $outputDir/dedup-report ]; then
						    mkdir $outputDir/dedup-report
						fi
	
						java \
						    -Xmx20g \
						    -Dsnappy.disable=true \
						    -jar /usr/local/bin/picard.jar MarkDuplicates \
						    REMOVE_DUPLICATES=true \
						    I=$outputDir/${strain}.${run}.sort.bam \
						    O=$outputDir/${strain}.${run}.sort.dedup.bam \
						    M=$outputDir/dedup-report/${strain}.${run}.txt \
						    VALIDATION_STRINGENCY=SILENT
						samtools index $outputDir/${strain}.${run}.sort.dedup.bam
						

					else
					    
					    echo "working on paired end data"
					    format=$(gzip -dc ${dir}/${strain}/${run}_pass_1.fastq.gz | head | /mnt/icy_1/inbredLines/software/DetermineFastqQualityEncoding.pl - 2> foo | grep "format" | grep -Eo "64|33")
					    echo "format = " $format
						if [ ${format} = "64" ]; then
							gzip -dc ${dir}/${strain}/${run}_pass_1.fastq.gz | /mnt/icy_1/inbredLines/software/convFQ.perl | gzip > ${dir}/${strain}/${run}_1.conv.fastq.gz
							rm 	${dir}/${strain}/${run}_pass_1.fastq.gz
							mv ${dir}/${strain}/${run}_1.conv.fastq.gz ${dir}/${strain}/${run}_1.fastq.gz
							
							gzip -dc ${dir}/${strain}/${run}_pass_2.fastq.gz | /mnt/icy_1/inbredLines/software/convFQ.perl | gzip > ${dir}/${strain}/${run}_2.conv.fastq.gz
							rm 	${dir}/${strain}/${run}_pass_2.fastq.gz
							mv ${dir}/${strain}/${run}_2.conv.fastq.gz ${dir}/${strain}/${run}_2.fastq.gz

						fi
					
						echo "mapping"
						bwa mem -t ${threads} \
						-R "@RG\tID:${strain};${run};pe\tSM:${strain}" \
						/mnt/internal_2/priscilla/ref_files/all_dmel.fasta \
						${dir}/${strain}/${run}_pass_1.fastq.gz \
						${dir}/${strain}/${run}_pass_2.fastq.gz | \
						samtools view -Suh - | \
						samtools sort -@ ${threads} -o ${outputDir}/${strain}.${run}.sort.bam
						
						samtools index ${outputDir}/${strain}.${run}.sort.bam

						echo "removing PCR duplicates"
						if [ ! -d $outputDir/dedup-report ]; then
						    mkdir $outputDir/dedup-report
						fi
	
						java \
						    -Xmx20g \
						    -Dsnappy.disable=true \
						    -jar /usr/local/bin/picard.jar MarkDuplicates \
						    REMOVE_DUPLICATES=true \
						    I=$outputDir/${strain}.${run}.sort.bam \
						    O=$outputDir/${strain}.${run}.sort.dedup.bam \
						    M=$outputDir/dedup-report/${strain}.${run}.txt \
						    VALIDATION_STRINGENCY=SILENT
						samtools index $outputDir/${strain}.${run}.sort.dedup.bam
						

					fi
					 
				done
							
			done
	
	
	}
	export -f mapReads_sraData

	parallel --gnu -j1 mapReads_sraData ::: "/mnt/internal_2/priscilla/HSparents/SRAdata"
