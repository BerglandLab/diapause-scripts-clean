
# Overview

The F4 and F5 hybrids were sequenced in three libraries, named PAE\_AOB\_01, PAE\_AOB\_02, and PAE\_AOB\_03. Within each library there are 960 barcoded samples, including some blank barcodes. The fastq files returned from the sequencing center are labeled with a number. The filenames.txt files correspond to the sample information in the order that they are numbered within each library (ie line 1 has the sample information for the files named *1\_R1\_001.fastq* and *1\_R2\_001.fastq*. The end result of these scripts are files *{sample\_name}\_merged.bam* and *{sample\_name}\_merged.bam.bai*.

## Map and clean up reads

* The reads were mapped with the pear\_bwa\_PAE\_AOB\_{1..3}.sh files
* PEAR combines paired ends into assembled and unassembled reads
* mapped with bwa and sorted with samtools
* PCR duplicates removed with Picard

# merge reads

* unassembled and assembled reads were merged with merge_bams.sh

