#thisscript.Rscript $draw $seed

library(data.table)
library(foreach)
library(doMC)
registerDoMC(16)
library(doRNG)

#make input file

args=commandArgs(trailingOnly=TRUE)

#draw and seed will be inputs to R script call
#draw is the imputation number (different for each set of randomly called heterozygotes)
draw=args[1]
print(draw)
#assign a  random seed in input filefor reproducibility

seed=as.numeric(args[2])
print(seed)
set.seed(seed)

#read in a file with "bad paths". These are haplotypes under 1 Mb that will be masked in final vcf. Has chr, start, stop, and haplotype (par1 or par2) that needs to be masked

bad.paths=fread("/scratch/pae3g/final_reconstruction2/bad_paths.csv")

#function to read in a single vcf and turn it into a diploid genotype wtih heterozygous parents randomly selected and short paths masked

import_data_mask<-function(fn){
    #pull out sample name from file name
  sample=strsplit(fn,split="[.]")[[1]][1]
  if (!file.size(fn) == 0) {
    #read in data and rename columns
    dat<-fread(paste('zcat ', fn, sep=""))
    print(sample)
    names(dat)=c("chr", "position", "h1", "h2")
    setkey(dat, chr, position)
    #pull out bad paths relevant to the file
    bad.paths.dat=bad.paths[sample.id==sample]
    #for each bad path segment, mask the appropriate haplotype by changing genotype to "./."
    foreach(bad.chr=bad.paths.dat$chr, bad.start=bad.paths.dat$cons.start, bad.stop=bad.paths.dat$cons.stop, bad.hap=bad.paths.dat$haplotype)%do%{
      if(bad.hap=="par1"){
        dat[chr==bad.chr & position>=bad.start & position<=bad.stop, h1:="./."]
      }
    if(bad.hap=="par2") {
      dat[chr==bad.chr & position>=bad.start & position<=bad.stop, h2:="./."]
      }
    }
    #randomly choose a homozygous parental genotype when parent is het
  dat[h1=="0/1", h1:=sample(c("0/0", "1/1"), 1)]
  dat[h2=="0/1", h2:=sample(c("0/0", "1/1"), 1)]
  #turn parental diploid genotypes into a single diploid genotype
  dat[h1=="0/0"&h2=="0/0",geno:="0|0"]
  dat[h1=="0/0"&h2=="1/1",geno:="0|1"]
  dat[h1=="1/1"&h2=="0/0",geno:="1|0"]
  dat[h1=="1/1"&h2=="1/1",geno:="1|1"]
  dat[h1=="0/0"&h2=="./.",geno:="0|."]
  dat[h1=="1/1"&h2=="./.",geno:="1|."]
  dat[h1=="./."&h2=="0/0",geno:=".|0"]
  dat[h1=="./."&h2=="1/1",geno:=".|1"]
  dat[h1=="./."&h2=="./.", geno:=".|."]
  #merge back to original vcf
  dat.universal<-merge(x, dat[,.(chr, position, geno)], all=TRUE)
  
  #make sure everything is in the right order
  dat.universal<-dat.universal[order(chr, position)]
  dat.universal[is.na(geno), geno:=".|."]
 
  setnames(dat.universal, "geno", sample)
  #print(dim(dat.universal))
  return(dat.universal[, -c("chr", "position", "ID", "REF", "ALT"), with=FALSE])
  }
  else {
    dat.universal=data.table(rep(".|.", dim(x)[1]))
    setnames(dat.universal, "V1", sample)
    return(dat.universal)
  }
}

#process samples for each swarm separately since they are in different folders
foreach(swarm=c("A", "B"))%do%{
  #set working directory for swarm
  wd=paste("/scratch/pae3g/genome-reconstruction/", swarm,  sep="")
  setwd(wd)
  file_list <- list.files(pattern = "PAE.*.14.vcf.gz$")
  #split files into 16 groups to prevent memory issues
  groups <- split(file_list, ceiling(seq_along(file_list)/ceiling(length(file_list)/16)))
  #read in the original vcf with all sites (just chr and pos)
  print("reading in master vcf...")
  x<-fread("/scratch/pae3g/genome-reconstruction/all_sites.vcf", header=TRUE)
  setnames(x, "#CHROM", "chr")
  setnames(x, "POS", "position")
  setkey(x, chr, position)
  #print(x)
  print("looping through grouped vcf files")

  foreach(i=(1:length(groups)))%dorng%{
      #use function to import each individual vcf and mask
    obj <- foreach(fn=groups[[i]], .combine="cbind")%do%{
      import_data_mask(fn)
    }
    #for first iteration, want to keep the initial vcf columns
    if (i==1 & swarm=="A"){
      obj=cbind(x, obj)
      setnames(obj, c("chr", 'position'), c("#CHROM", "POS"))
      obj[,QUAL:=rep(1, dim(obj)[1])]
      obj[,FILTER:=rep(NA, dim(obj)[1])]
      #save seed in info column
      obj[,INFO:=rep(seed, dim(obj)[1])]
      obj[,FORMAT:=rep("GT", dim(obj)[1])]
      #order in traditional vcf order
      setcolorder(obj, c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", names(obj)[6:(length(obj)-4)]))

      write.table(obj, paste("/scratch/pae3g/genome-reconstruction/", swarm,"het_random_draw", draw,  "_part", i, ".vcf", sep=""), sep="\t", row.names=FALSE, quote=FALSE)
    } else {
        #don't add all the extra vcf columns for parts 2-16
      write.table(obj, paste("/scratch/pae3g/genome-reconstruction/", swarm,"het_random_draw", draw,  "_part", i, ".vcf", sep=""), sep="\t", row.names=FALSE, quote=FALSE)
    }
    rm(obj)
    gc()
  }
}


