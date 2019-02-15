##prepare files
prepare_read_simulation_fasta.sh takes founders vcf file, splits it into individual vcfs, then makes fasta files for each chromosome of each individual, using ambiguous bases for heterozygous sites


##generate simulated hybrid swarm
nannyGenerateMapFreqs runs the GenerateMapFreqs.slurm script for different parameters to generate simulated F4 and F5 hybrid swarms from the founding lines

##reconstruct simulated genomes
extract_simulated_bams.sh extracts the bam files produced for genome reconstruction

reconstruct_sim_PAE.14.sh runs the reconstruction on simulated bams

#calculate accuracy

calc_accuracy_array.sh compares the genotypes of reconstructed genomes to the true sequence that went into the read simulator

