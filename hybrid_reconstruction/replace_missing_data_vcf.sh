#!/bin/bash

draw=${1}

#first used R to randomly choose an allele if parent is heterozygous

#next stitch these files together as separate swarms


    echo "working on draw " $draw
    for swarm in A B; do
	echo "combining vcf parts swarm " $swarm
	if [[ ! -f  /scratch/pae3g/genome-reconstruction/$swarm"_het_random_draw"$draw.vcf ]]; then
	    paste -d "\t" /scratch/pae3g/genome-reconstruction/$swarm"het_random_draw"$draw"_part"{1..16}.vcf > /scratch/pae3g/genome-reconstruction/$swarm"_het_random_draw"$draw.vcf
	fi
	

	#run awk script to sub in missing data on two files separately (so that most common allele for each population is chosen) and make separate logs for each that can be used for later filtering

	echo "awk replacing missing data " $swarm
	awk -v swarm=$swarm -v draw=$draw -v OFS="\t" 'BEGIN {
    print "chr\tpos\tnRef\tnAlt\tmostLikelyGeno\tn\ttotal">"/scratch/pae3g/genome-reconstruction/replace_"swarm"_"draw"_better.log"
}
{
    nRef=0
    nAlt=0
    for(i=1; i<=NF; i++) {
	if($i=="0|0") nRef+=2 
	if($i=="0|." || $i==".|0") nRef+=1
	if($i=="1|." || $i==".|1") nAlt+=1
	if($i=="0|1" || $i=="1|0"){
	    nRef+=1
	    nAlt+=1
	}
	if($i=="1|1") nAlt+=2
    }
    if((nRef+nAlt)==0){
	refFreq=1
	altFreq=0
    }
    else{
	refFreq=nRef/(nRef+nAlt)
	altFreq=1-refFreq
    }
    freq00=refFreq*refFreq
    freq01=refFreq*altFreq*2
    freq11=altFreq*altFreq
    if(freq00 > freq01 && freq00 > freq11){
	mostLikelyGeno="homRef"
	n=gsub(/\./, 0); print
	}
    if(freq01 > freq00 && freq01 > freq11){
	mostLikelyGeno="het"
	n=gsub(/0\|\./, "0|1") + gsub(/\.\|0/, "1|0") + gsub(/1\|\./, "1|0") + gsub(/\.\|1/, "0|1")+gsub(/\.\|\./, "0|1")*2; print
    }
    if(freq11 > freq01 && freq11>freq00){
	mostLikelyGeno="homAlt"
	n=gsub(/\./, 1); print
    }
total=nRef+nAlt+n
{  printf $1"\t"$2"\t"nRef"\t"nAlt"\t"mostLikelyGeno"\t"n"\t"total"\n" > "/scratch/pae3g/genome-reconstruction/replace_"swarm"_"draw"_better.log"}
}' /scratch/pae3g/genome-reconstruction/$swarm"_het_random_draw"$draw.vcf > /scratch/pae3g/genome-reconstruction/$swarm"_het_random_draw"$draw"_replaced.vcf"

    done

    echo "combining indivdual swarm vcfs "
    paste /scratch/pae3g/genome-reconstruction/A_het_random_draw$draw"_replaced.vcf" /scratch/pae3g/genome-reconstruction/B_het_random_draw$draw"_replaced.vcf" > /scratch/pae3g/genome-reconstruction/final2_draw$draw"_replaced.vcf"


