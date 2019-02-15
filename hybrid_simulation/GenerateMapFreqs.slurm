#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=5
#SBATCH --mem 30000
#SBATCH -t 0-01:30:00
#SBATCH -p standard 
#SBATCH --account bergland-erickson

# Copy and paste when running manually:

# module load R
# topDirectory="/scratch/$USER/genome-reconstruction/"
# referenceGroup=dgrp2
# referenceGenome="$topDirectory"/variants/$referenceGroup/ref/all_dmel.fasta
# nGenerations="5"
# nFounders="32"
# populationSize="100"
# growthRate="1"
# bedfile=/"$topDirectory"/variants/$referenceGroup/recombination.bed
# sex="dioecious"
# recombination="femaleOnly"
# nThreads="1"
# vcf="$topDirectory"/variants/$referenceGroup/vcf/dgrp2.2L.vcf
# nIndsPerPop="2"
# readLength=100
# coverage=0.05
# rep=234
# ind=77

module load R/3.3.0


topDirectory=$1
referenceGroup=$2
referenceGenome=$3
nGenerations=$4
nFounders=$5
populationSize=$6
growthRate=$7
bedfile=$8
sex=$9
recombination=${10}
nThreads=${11}
vcf=${12}
nIndsPerPop=${13}
readLength=${14}
coverage=${15}
rep=${16}

echo 1 $1
echo 2 $2
echo 3 $3
echo 4 $4
echo 5 $5
echo 6 $6
echo 7 $7
echo 8 $8
echo 9 $9
echo 10 ${10}
echo 11 ${11}
echo 12 ${12}
echo 13 ${13}
echo 14 ${14}
echo 15 ${15}
echo 16 ${16}
## FUNCTIONS

function getNReads {

# getNReads $readLength $coverage $genomeSize

Rscript - <<-EOF
    readLength <- as.numeric($1)
    coverage <- as.numeric($2)
    genomeSize <- as.numeric($3)
    nReads = round(genomeSize * coverage / (4*readLength))
    cat(nReads)
EOF
}

export -f getNReads

function subsetHarpPriorsCSV {

priorsCSV=$1
foundersFile=$2

Rscript - <<-EOF
    
    library(data.table)
    
    priors.file <- "$priorsCSV"
    founders.file <- "$foundersFile"
    
    priors <- fread(priors.file, sep=",", header=TRUE, showProgress=FALSE)
    founders <- readLines(founders.file)
    
    columnsToKeep <- c(colnames(priors)[1:2], founders, "Coverage")
    write.table(priors[,columnsToKeep, with=FALSE], file="", quote=FALSE, col.names=TRUE, row.names=FALSE, sep=",")
EOF
}

export -f subsetHarpPriorsCSV

function extractVCF {
python - <<-EOF

sites_filename = "$topDirectory/variants/$referenceGroup/vcf/$referenceGroup.$chr.vcf.sites"

with open(sites_filename, 'r') as infile:
    sites = infile.readlines()

sites = [int(x.rstrip()) for x in sites]

def getHaplotype(fasta_filename, sites):
    with open(fasta_filename, 'r') as infile:
        fasta = infile.read().rstrip()
        return([fasta[x-1] for x in sites])

Haplotype1 = getHaplotype("Chr$chr.1.fasta.tmp", sites)
Haplotype2 = getHaplotype("Chr$chr.2.fasta.tmp", sites)

print "\n".join(["$chr\t"+"\t".join(x) for x in zip([str(x) for x in sites], Haplotype1, Haplotype2)])

EOF
}

export -f extractVCF

function GenerateMapFreq {

    topDirectory=$1
    tmpWorkDir=$2
    referenceGenome=$3
    readLength=$4
    coverage=$5
    iteration=$6
    ind=$7
    referenceGroup=$8

    # WGSIM PIPELINE

    wgsim="$topDirectory"/etc/wgsim-master/wgsim
    picard="$topDirectory"/etc/picard.jar
    pear="$topDirectory"/etc/pear
    bwa="$topDirectory"/etc/bwa/bwa
    samtools="$topDirectory"/etc/samtools
    harp="$topDirectory"/etc/harp
    gatk="$topDirectory"/etc/GenomeAnalysisTK.jar
    jobDir="$tmpWorkDir"/"$iteration"_"$ind"/
    window_step=100000
    window_width=100000
    
    # Load chromosome lengths from saved associative array ${lengths[@]}
    echo loading chromosome length
    source "$topDirectory"/variants/$referenceGroup/seq/lengths.sh

    mkdir -p $jobDir
    cp "$topDirectory"/output/$referenceGroup/$iteration/population.haps $jobDir
    cp "$topDirectory"/output/$referenceGroup/$iteration/founders.txt $jobDir
    cd $jobDir
    
    while read chromosome haplotype start stop indNo sex lineID; do
        # Adjust cut command to correctly point to .seq file names based on .haps file lineID
        printf '%s' $(cut -c"$start"-"$stop" "$topDirectory"/variants/$referenceGroup/seq/"$chromosome"/"$lineID"_Chr"$chromosome".seq) >> Chr$chromosome.$haplotype.fasta.tmp
    done < <(awk -F"\t" -v indNo=$ind 'NR >1 && $5 == indNo' population.haps)
    
    # if male; duplicate X haplotype to second file
    if [ ! -f "ChrX.1.fasta.tmp" ]; then
        cp ChrX.2.fasta.tmp ChrX.1.fasta.tmp
    fi
    
    if [ ! -f "ChrX.2.fasta.tmp" ]; then
        cp ChrX.1.fasta.tmp ChrX.2.fasta.tmp
    fi
    
    for chr in ${!lengths[@]}; do
        extractVCF >> $ind.vcf.reduced
    done
    
    
    # Add header lines
    for file in *.fasta.tmp; do
        # replace 101+ repeats of N with exactly 100 Ns
        #sed -i 's/N\{501,\}/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/g' $file
        # add headers
            sed -i "1 i\>${file%.fasta}" $file
        
        # add newline to end of file
            sed -i -e '$a\' $file
    done
    
    cat Chr*.fasta.tmp | fold > sim.fasta && rm *.tmp

    
    genomeSize=$(grep "^[^>]" sim.fasta | wc -c)

nReads=$(getNReads $readLength $coverage $genomeSize)

    # Generate fastq files with wgsim
        echo generating reads with wgsim
        $wgsim -1 $readLength -2 $readLength -N $nReads -e 0.001 -r 0 -R 0 sim.fasta $ind.F.fq $ind.R.fq && rm sim.fasta && \
        bzip2 *.fq
        
# MAP READS PIPELINE

    $pear \
    -f $ind.F.fq.bz2 \
    -r $ind.R.fq.bz2 \
    -o $ind
    
    echo 'mapping assembled reads with bwa and exporting to samtools'
    $bwa mem -t 1 \
            -R "@RG\tID:simID;simRun;simLane;assembled\tSM:simSM" \
            $referenceGenome \
            $ind.assembled.fastq | \
            $samtools view -Suh - | \
            $samtools sort -o $ind.assembled.sort.bam
        echo "samtools sorting completed, now indexing"
            $samtools index $ind.assembled.sort.bam
            echo "removing PCR duplicates"

            java \
            -Xmx2G \
            -Dsnappy.disable=true \
            -jar $picard  MarkDuplicates \
            REMOVE_DUPLICATES=true \
            I=$ind.assembled.sort.bam \
            O=$ind.assembled.sort.dedup.bam \
            M=$ind.dups
            $samtools index $ind.assembled.sort.dedup.bam

    echo 'mapping unassembled reads with bwa and exporting to samtools'
    $bwa mem -t 1 \
            -R "@RG\tID:simID;simRun;simLane;unassembled\tSM:simSM" \
            $referenceGenome \
            $ind.unassembled.forward.fastq \
            $ind.unassembled.reverse.fastq | \
            $samtools view -Suh - | \
            $samtools sort -o $ind.unassembled.sort.bam
        echo "samtools sorting completed, now indexing"
            $samtools index $ind.unassembled.sort.bam
            echo "removing PCR duplicates"

            java \
            -Xmx2G \
            -Dsnappy.disable=true \
            -jar $picard  MarkDuplicates \
            REMOVE_DUPLICATES=true \
            I=$ind.unassembled.sort.bam \
            O=$ind.unassembled.sort.dedup.bam \
            M=$ind.dups
            $samtools index $ind.unassembled.sort.dedup.bam

        $samtools merge $ind.bam $ind.assembled.sort.dedup.bam $ind.unassembled.sort.dedup.bam && \
        $samtools index $ind.bam && \
        rm -rf $ind.tmp && \
        rm $ind.F.fq.bz2 && \
        rm $ind.R.fq.bz2 && \
        rm $ind.dups && \
        rm $ind.discarded.fastq && \
        rm *assembled*
        
        # Record number of reads mapped per chromsome
        $samtools idxstats $ind.bam > $ind.bam.stats

# GET READ COUNTS HERE




# READ COUNTS & HARP FREQ PIPELINE

    # iterate through chromosomes; 
    for chromosome in ${!lengths[@]}; do
    
        java \
        -Xmx2G \
        -jar $gatk \
        -R $referenceGenome \
        -T ASEReadCounter \
        -I $ind.bam \
        -o $ind.$chromosome.readcounts \
        -sites $topDirectory/variants/$referenceGroup/vcf/$referenceGroup.$chromosome.vcf
        
        echo running harp on chromosome $chromosome
        
        mkdir -p $jobDir/$chromosome
        cd $jobDir/$chromosome
        length=${lengths[$chromosome]}
        
        # Subset priors, writing to priors.csv for this specific chromosome folder
        subsetHarpPriorsCSV \
            $topDirectory/variants/$referenceGroup/priors/prior_$chromosome.noRep.noIndel.noHet.csv \
            $topDirectory/output/$referenceGroup/$iteration/founders.txt \
            > priors.csv
        
        # index harp snp file
        $topDirectory/etc/index_snp_table priors.csv 1
        
        # Run harp like
        echo running harp like
        $harp like \
        --bam $jobDir/$ind.bam \
        --region $chromosome:1-$length \
        --refseq $referenceGenome \
        --snps priors.csv \
        --stem $ind.$chromosome
    
        # Run harp freq
        echo running harp freq
        $harp freq \
        --bam $jobDir/$ind.bam \
        --region $chromosome:1-$length \
        --refseq $referenceGenome \
        --snps priors.csv \
        --stem $ind.$chromosome \
        --window_step $window_step \
        --window_width $window_width \
        --em_min_freq_cutoff 0.0001
    
        # Cleanup 
        echo finished with "$chromosome" freqs, cleaning up and saving output
        mv "$ind."$chromosome.freqs $jobDir && \
        cd $jobDir && \
        rm -rf $jobDir/$chromosome
    done
    
    # gzip everything to saved output > iteration.tar.gz
    cd $jobDir && \
    rm population.haps && \
    rm founders.txt && \
    tar -cvf - . | gzip > $topDirectory/output/$referenceGroup/$iteration/$ind.tar.gz
    cd $topDirectory && \
    rm -rf $jobDir
}

export -f GenerateMapFreq


## generate forward-simulated population

cd "$topDirectory"/output/
iteration="$nFounders"F_"$nGenerations"G_"$coverage"X_"$populationSize"N_rep"$rep"

# change into directory for output

mkdir -p "$topDirectory"/output/$referenceGroup/$iteration 
cd "$topDirectory"/output/$referenceGroup/$iteration
tmpWorkDir=/dev/shm/$USER/$JOBID/


Rscript "$topDirectory"/scripts/forwardSimulator.Rscript \
"-bed" "$bedfile" \
"-nfounders" "$nFounders" \
"-ngenerations" "$nGenerations" \
"-n0" "$populationSize" \
"-sex" "$sex" \
"-rate" "$growthRate" \
"-recombination" "$recombination" \
"-vcf" "$vcf" \
"-nthreads" "$nThreads"


parallel -j $nThreads GenerateMapFreq ::: $topDirectory ::: $tmpWorkDir ::: $referenceGenome ::: $readLength ::: $coverage ::: $iteration ::: $(seq 1 $nIndsPerPop) ::: $referenceGroup

