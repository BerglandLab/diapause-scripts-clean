#!/bin/bash


parameterFile="/scratch/$USER/genome-reconstruction/parameters_sim.txt"

reconstructionGroup=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $parameterFile | cut -f 1 )
ind_id=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $parameterFile | cut -f 2 )
nGenerations=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $parameterFile | cut -f 3 )

module load R/3.3.0
cd ${workDir}

Rscript - ${reconstructionGroup} ${ind_id} ${nGenerations} <<EOF > /scratch/pae3g/genome-reconstruction/accuracy/$ind_id.summary.dat
#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(foreach)

args <- commandArgs(trailingOnly=TRUE)
reconstructionGroup <- args[1]
ind_id <- args[2]
nGenerations <- args[3]


#for testing
#ind_id="34F_4G_0.5X_10000N_rep100.1.A"
#nGenerations=4
#reconstructionGroup="A"


estimate.filename <- paste("zcat  /scratch/pae3g/genome-reconstruction/", reconstructionGroup, "/", ind_id, ".estimate.vcf.gz", sep="")
estimate.filename.14 <- paste("zcat /scratch/pae3g/genome-reconstruction/", reconstructionGroup, "/", ind_id, ".estimate.14.vcf.gz", sep="")

true.vcf.filename <- paste("/scratch/pae3g/genome-reconstruction/output/cage", reconstructionGroup,"/", ind_id, ".vcf.reduced", sep="")
sites.filename <- "/scratch/pae3g/genome-reconstruction/sites.vcf"

estimate.vcf <- fread(estimate.filename, header=FALSE, na.strings=c("./.", "0/1"), showProgress=FALSE, col.names=c("CHR","POS","e1","e2"), key=c("CHR","POS"))
estimate.vcf.14 <- fread(estimate.filename.14, header=FALSE, na.strings=c("./.","0/1"), showProgress=FALSE, col.names=c("CHR","POS","e1.14","e2.14"), key=c("CHR","POS"))

true.vcf <- fread(true.vcf.filename, header=TRUE, na.strings="./.", showProgress=FALSE, col.names=c("CHR", "POS","a1","a2"))
sites.vcf    <- fread(sites.filename, header=TRUE, na.strings="./.", showProgress=FALSE, col.names=c("CHR","POS", "ref", "alt"), key=c("CHR","POS"))

# combine true chromosome/site columns with true genotype columns
true.vcf <- merge(sites.vcf, true.vcf, by=c("CHR", "POS"))

# combine truth with estimate
all.vcf <- merge(true.vcf, estimate.vcf, by=c("CHR", "POS"))
all.vcf <- merge(all.vcf, estimate.vcf.14, by=c("CHR", "POS"))

# remove any rows containing missing information
all.vcf <- all.vcf[!is.na(a1) & ! is.na(a2) & ! is.na(e1) & ! is.na(e2) & ! is.na(e1.14) & ! is.na(e2.14)]

#convert actual calls to 0/1
all.vcf[,a1:=ifelse(a1==ref, 0,1)]
all.vcf[,a2:=ifelse(a2==ref, 0,1)]


# Convert "0/0" and "1/1" to "0" and "1"
if(any(all.vcf=="0/0") | any(all.vcf=="1/1")) {
  all.vcf[, c("e1", "e2", "e1.14", "e2.14") := lapply(.SD, factor, levels=c("0/0","1/1")), .SDcols=c("e1", "e2", "e1.14", "e2.14")]
  all.vcf[, c("e1", "e2", "e1.14", "e2.14") := lapply(.SD, as.numeric), .SDcols=c("e1", "e2", "e1.14", "e2.14")]
  all.vcf[, c("e1", "e2", "e1.14", "e2.14") := lapply(.SD, "-", 1), .SDcols=c("e1", "e2", "e1.14", "e2.14")]
}

# calculate alt dosage by summing both "actual" alleles, and summing both "estimate" alleles
# "actual" and "estimate" columns will take on values of either (0, 1, 2) where 0=ref/ref, 1=heterozygote, 2=alt/alt
all.vcf[, actual := apply(.SD, 1, function(x) sum(x == 1)), .SDcols=c("a1","a2")]
all.vcf[, estimate := apply(.SD, 1, function(x) sum(x == 1)), .SDcols=c("e1","e2")]
all.vcf[, estimate.14 := apply(.SD, 1, function(x) sum(x == 1)), .SDcols=c("e1.14","e2.14")]

# Iterate through chromosomes, calculating  the fraction of sites where actual == estimate
dat.all <- rbindlist(foreach(chr.i = unique(all.vcf[,CHR])) %do% {
  nSitesVariableInPopulation <- dim(sites.vcf[.(chr.i)])[1]
  nSites <- dim(all.vcf[.(chr.i)])[1]
  nMatches <- dim(all.vcf[.(chr.i)][actual==estimate])[1]
  percentMatches <- nMatches/nSites
  nMatches.14 <- dim(all.vcf[.(chr.i)][actual==estimate.14])[1]
  percentMatches.14 <- (nMatches.14)/nSites
  return(data.table("CHR"=chr.i, nSitesVariableInPopulation, nSites, nMatches, percentMatches, nMatches.14, percentMatches.14))
})

# add parameter columns
dat.all[, "ind_id" := ind_id]
dat.all[, "nGenerations" := nGenerations]

# write output
setcolorder(dat.all, c("ind_id","nGenerations", "CHR", "nSitesVariableInPopulation", "nSites", "nMatches", "percentMatches", "nMatches.14", "percentMatches.14"))
write.table(dat.all, file=paste("/scratch/pae3g/genome-reconstruction/accuracy/", ind_id, ".summary.dat", sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

EOF
