####Parental Strain Processing

###parental sequencing data can be found in NCBI Bioproject PRJNA522357


###download sra sequences in sraData.delim

### map new sequence reads with pear_bwa.sh

###map SRA reads with map_sra2.sh


###one DGRP file is especially large adn causes problems in
GATK. downsample it
samtools view -s .2 -b DGRP_28243.SRR835347.sort.dedup.renamed.bam > DGRP_28243.SRR835347.sort.dedup.renamed.downsampled.bam

###merge together new bams and old bams, except for 3 samples where
the new sequence does not match the old sequence

###use haplotype_caller_rivanna.sh to make gVCFs for each of these
merged samples





