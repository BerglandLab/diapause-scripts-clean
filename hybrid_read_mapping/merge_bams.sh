#!/bin/bash
#thisscript bam_directory file_list.txt

mergeBams () {
directory=$1
file=$2


    afile=$directory/$file.assembled.sort.dedup.bam
    echo $filename
    ufile=$directory/$file.unassembled.sort.dedup.bam
    samtools merge $directory/$file.merged.bam $afile $ufile
    samtools index $directory/$file.merged.bam 

}

export -f mergeBams

files=($(cat /mnt/pricey_2/priscilla/hybrid/raw_data/PAE_AOB_01/PAE_AOB_01_filenames.txt))

parallel --gnu -j10 mergeBams ::: "/mnt/pricey_2/priscilla/hybrid/mapped/PAE_AOB_01" ::: "${files[@]}"

files=($(cat /mnt/pricey_2/priscilla/hybrid/raw_data/PAE_AOB_02/PAE_AOB_02_filenames.txt))

parallel --gnu -j10 mergeBams ::: "/mnt/pricey_2/priscilla/hybrid/mapped/PAE_AOB_02" ::: "${files[@]}"

files=($(cat /mnt/pricey_2/priscilla/hybrid/raw_data/PAE_AOB_03/PAE_AOB_03_filenames.txt))

parallel --gnu -j10 mergeBams ::: "/mnt/pricey_2/priscilla/hybrid/mapped/PAE_AOB_03" ::: "${files[@]}" 
