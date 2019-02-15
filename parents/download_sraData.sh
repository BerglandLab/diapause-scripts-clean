#!/bin/bash

downloadSRA () {
	echo ${1} 
	
	srr=$( echo ${1} | cut -f1 -d'\t' )
	#outDir=/mnt/spicy_2/dest/00_fastq
	fileStem=$( echo {1} | cut -f2, -d '\t' )
	
	firstThree=${srr:0:6}

	#pop=$( echo $outDir | cut -f1 -d'/' ) 

	#if [ ! -d /mnt/inbred/inbredLines/00_rawData/${pop} ]; then
	#	mkdir /mnt/inbred/inbredLines/00_rawData/${pop}
	#fi
	#
	#if [ ! -d /mnt/inbred/inbredLines/00_rawData/${outDir} ]; then
	#	mkdir /mnt/inbred/inbredLines/00_rawData/${outDir}
	#fi

	/home/bergland/.aspera/connect/bin/ascp -QTr -l 400M \
	-i /home/bergland/.aspera/connect/etc/asperaweb_id_dsa.openssh \
	anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/${firstThree}/${srr}/${srr}.sra \
	/mnt/spicy_2/dest/00_fastq/${fileStem}.sra

	fastq-dump --split-files \
	--gzip \
	-O /mnt/inbred/inbredLines/00_rawData/ \
	/mnt/spicy_2/dest/00_fastq/${fileStem}.sra

	#rm /mnt/inbred/inbredLines/00_rawData/${outDir}/${srr}.sra

}
export -f downloadSRA

parallel --gnu -j1 -a /mnt/inbred/inbredLines/scripts/sraData.delim downloadSRA 


























## SE/Carribean
#	wget -O ./PRJNA274815_info.csv 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term= PRJNA274815'
#	cat PRJNA274815_info.csv | cut -f1,14 -d',' | sed 's/Dmel WGS //g' | sed 's/ Experiment /,/g' | awk -F',' '
#	BEGIN {
#		id[1]="20_28"
#		id[2]="20_17"
#		id[3]="13_34"
#		id[4]="13_29"
#		id[5]="4_12"
#		id[6]="4_27"
#		id[7]="21_39"
#		id[8]="21_36"
#		id[9]="24_2"
#		id[10]="24_9"
#		id[11]="28_8"
#		id[12]="33_16"
#		id[13]="33_11"
#		id[14]="36_9"
#		id[15]="36_12"
#		id[16]="40_23"
#		id[17]="40_10"
#		id[18]="42_23"
#		id[19]="42_20"
#		id[20]="43_19"
#		id[21]="43_18"
#		id[22]="H_29"
#		id[23]="H_25"
#		
#		newName[1]="Selba_AL_20_28"
#		newName[2]="Selba_AL_20_17"
#		newName[3]="Thomasville_GA_13_34"
#		newName[4]="Thomasville_GA_13_29"
#		newName[5]="TampaBay_FL_4_12"
#		newName[6]="TampaBay_FL_4_27"
#		newName[7]="Birmingham_AL_21_39"
#		newName[8]="Birmingham_AL_21_36"
#		newName[9]="Meridian_MS_24_2"
#		newName[10]="Meridian_MS_24_9"
#		newName[11]="Sebastian_FL_28_8"
#		newName[12]="Freeport_GrandBahamasWest_33_16"
#		newName[13]="Freeport_GrandBahamasWest_33_11"
#		newName[14]="GeorgeTown_Exumas_36_9"
#		newName[15]="GeorgeTown_Exumas_36_12"
#		newName[16]="BullocksHarbor_BerryIslands_40_23"
#		newName[17]="BullocksHarbor_BerryIslands_40_10"
#		newName[18]="Cockburn_SanSalvador_42_23"
#		newName[19]="Cockburn_SanSalvador_42_20"
#		newName[20]="Mayaguana_Mayaguana_43_19"
#		newName[21]="Mayaguana_Mayaguana_43_18"
#		newName[22]="PortAuPrince_Haiti_H_29"
#		newName[23]="PortAuPrince_Haiti_H_25"
#	}	
#		
#	{
#	
#		for(i=1; i<=23; i++) {
#			if(id[i]==$2) {
#				print ""newName[i]"/"$3"/ "$1
#			}
#		}
#	}'
#		
	
		
	
