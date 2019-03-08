#Overview

These scripts process new sequencing data for the inbred lines used
to initiate the hybrid swarms. The new sequencing data was combined
with existing sequencing data (available on the SRA) for variant
calling. Two vcfs are produced: one containing the most
high-confidence variants to be used for reconstruction, and a second
vcf with more variants used for association mapping.

## Download and process existing parental sequence reads

Use sraData.delim and download sraData.sh to download appropriate
sequence files. Use mapr_sra2.sh to map and clean up the
reads. Because one DGRP file is particularly large, the mapped reads
were downsampled with following command: `samtools view -s .2 -b
DGRP_28243.SRR835347.sort.dedup.renamed.bam >
DGRP_28243.SRR835347.sort.dedup.renamed.downsampled.bam`


## merge existing data and new data bam files


###use haplotype_caller_rivanna.sh to make gVCFs for each of these
merged samples





