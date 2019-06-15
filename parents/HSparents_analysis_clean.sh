
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
