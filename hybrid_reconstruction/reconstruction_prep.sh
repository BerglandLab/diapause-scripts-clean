#move files to rivanna cluster

## filter parental snps for reconstruction using files produced in filter_parents.R

bedtools intersect -header -a /mnt/pricey_2/priscilla/hybrid/hc/hs.hc.A.vcf -b /mnt/pricey_2/priscilla/hybrid/hc/hs.hc.A.keep.bed > /mnt/pricey_2/priscilla/hybrid/hc/hs.hc.A.keep.vcf

bedtools intersect -header -a /mnt/pricey_2/priscilla/hybrid/hc/hs.hc.B.vcf -b /mnt/pricey_2/priscilla/hybrid/hc/hs.hc.B.keep.bed > /mnt/pricey_2/priscilla/hybrid/hc/hs.hc.B.keep.vcf


#haplotypes.vcf has the highest quality snps for genome reconstruction
scp hs.hc.A.keep.vcf pae3g@rivanna1.hpc.virginia.edu:/scratch/pae3g/genome-reconstruction/A/haplotypes.vcf
scp hs.hc.B.keep.vcf pae3g@rivanna1.hpc.virginia.edu:/scratch/pae3g/genome-reconstruction/B/haplotypes.vcf

#reconstruct.vcf has all genotypes to prepare final hybrid vcfs
scp hs.hc.snp.99.9.A.vcf pae3g@rivanna1.hpc.virginia.edu:/scratch/pae3g/genome-reconstruction/A/reconstruct.vcf
scp hs.hc.snp.99.9.B.vcf pae3g@rivanna1.hpc.virginia.edu:/scratch/pae3g/genome-reconstruction/B/reconstruct.vcf

#on rivanna, make bgzip and .tbi files in A and B folders:

/scratch/$USER/genome-reconstruction/etc/htslib/bin/bgzip -c haplotypes.vcf > haplotypes.vcf.gz && /scratch/$USER/genome-reconstruction/etc/htslib/bin/tabix -p vcf haplotypes.vcf.gz

cp /scratch/pae3g/hybrid/seqs/all_dmel* /scratch/pae3g/genome-reconstruction/A
cp /scratch/pae3g/hybrid/seqs/all_dmel* /scratch/pae3g/genome-reconstruction/B

    # Create a file "founders.txt" within this directory, containing a list of founding lines.
zgrep -m 1 "^#CHROM" haplotypes.vcf.gz | cut -f 10- | tr "\t" "\n" > founders.txt

#make priors file for harp
function subsetPriors {
    chromosome=${1}
    priorHeader=$(zcat haplotypes.vcf.gz | grep -m 1 "^#CHROM" | cut -f 10- | tr "\t" ","),Coverage
    echo $chromosome,Ref,$priorHeader, > $chromosome.priors.csv

    zgrep -w "^$chromosome" haplotypes.vcf.gz | python -c '
import sys
for line in sys.stdin:
    if line.startswith("#CHROM"):
        splitLine = line.rstrip().split()
        ids = splitLine[9:]
        output = ["$chromosome", "Ref"] +  ids + ["Coverage,"]
        print ",".join(output)
    elif not line.startswith("##"):
        splitLine = line.rstrip().split()
        chromosome, position, id, ref, alt = splitLine[0:5]
        output = [position, ref] + [ref if x in ("0/0", "0") else alt if x in ("1/1", "1")  else "N" for x in splitLine[9:]]
        print ",".join(output)+",34,"' >> $chromosome.priors.csv
}
export -f subsetPriors

chromosomes=( 2L 2R 3L 3R X )
parallel -j 1 subsetPriors ::: ${chromosomes[@]}

#make reconstruction.vcf.gz for A and B

/scratch/$USER/genome-reconstruction/etc/htslib/bin/bgzip -c reconstruct.vcf > reconstruct.vcf.gz && /scratch/$USER/genome-reconstruction/etc/htslib/bin/tabix -p vcf reconstruct.vcf.gz

cd ../A
/scratch/$USER/genome-reconstruction/etc/htslib/bin/bgzip -c reconstruct.vcf > reconstruct.vcf.gz && /scratch/$USER/genome-reconstruction/etc/htslib/bin/tabix -p vcf reconstruct.vcf.gz
