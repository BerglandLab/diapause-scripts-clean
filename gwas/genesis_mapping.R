#!/usr/bin/env Rscript

library(GWASTools, quietly=TRUE)
library(GENESIS, quietly=TRUE)
library(data.table, quietly=TRUE)
library(SNPRelate, quietly=TRUE)
library(gdsfmt, quietly=TRUE)
library(doMC, quietly=TRUE)
registerDoMC(20)
library(foreach, quietly=TRUE)

#make input data for this analysis
#input<-data.table(GRM=rep(c("Eigenstrat", "Corr", "EIGMIX", "GCTA", "IndivBeta", "KING-robust"), each=9), pop=rep(c("A", "B", "both"), each=1, times=18), phenotype=rep(c("diapause.bin", "diapause.bin9", "eggP"), each=3, times=6))
#write.table(input, "/scratch/pae3g/final_reconstruction2/genesis_input.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

args=commandArgs(trailingOnly=TRUE)
filestem=args[1]
pop=args[2]
phenotype=args[3]
draw=args[4]
perm=args[5]
seed=args[6]


### modified functions
GdsGenotypeReader.par <- function (filename, genotypeDim, genotypeVar, snpIDvar, scanIDvar, allow.fork=T, ...) {
  if (missing(filename)) 
    stop("filename is required")
  if (missing(genotypeVar)) 
    genotypeVar <- "genotype"
  if (missing(snpIDvar)) 
    snpIDvar <- "snp.id"
  if (missing(scanIDvar)) 
    scanIDvar <- "sample.id"
  input.gds <- is(filename, "gds.class")
  tmpobj <- GdsReader.par(filename, allow.fork=allow.fork)
  snpDim <- getDimension(tmpobj, snpIDvar)
  scanDim <- getDimension(tmpobj, scanIDvar)
  genoDim <- getDimension(tmpobj, genotypeVar)
  if (missing(genotypeDim)) {
    if (snpDim == scanDim) {
      genotypeDim <- ""
    }
    else if (all(genoDim == c(snpDim, scanDim))) {
      genotypeDim <- "snp,scan"
    }
    else if (all(genoDim == c(scanDim, snpDim))) {
      genotypeDim <- "scan,snp"
    }
    else {
      genotypeDim <- ""
    }
  }
  tryCatch(new("GdsGenotypeReader", tmpobj, genotypeDim = genotypeDim, 
               genotypeVar = genotypeVar, snpIDvar = snpIDvar, scanIDvar = scanIDvar, 
               ...), error = function(e) {
                 if (!input.gds) 
                   close(tmpobj)
                 stop(e)
               })
}

GdsReader.par <- function (filename, allow.fork) {
  if (missing(filename)) 
    stop("filename is required")
  if (is(filename, "gds.class")) {
    input.gds <- TRUE
    handler <- filename
    filename <- handler$filename
  }
  else {
    input.gds <- FALSE
    if (!file.exists(filename)) 
      stop("Error in opening file ", filename, ": no such file or directory")
    handler <- openfn.gds(filename, allow.fork=allow.fork)
  }
  new("GdsReader", filename = filename, handler = handler)
}


### genoData 	
geno <- GdsGenotypeReader.par(filename = paste(filestem, ".vcf.gds", sep=""), allow.fork=T)
genoData <- GenotypeData(geno)

  #read data
phenos <- fread("/scratch/pae3g/phenos/phenos_062018.txt")
setkey(phenos, id)
  ### some twiddles	
  phenos[,sex := "F"]	
  if (pop=="both"){
    load(paste(filestem, "_LOCO_Eigenstrat_strict_both.Rdat", sep=""))
    } else {
    load(paste("/scratch/pae3g/final_reconstruction2/", GRM, "_LOCO_GRM_replaced", pop, ".Rdat", sep=""))
  }
  
  phenos <- phenos[id %in% rownames(loco.list[[1]])]
  #if a permutation, randomly scramble sample ids within a swarm. if perm is 0, ids will stay the same
  if(perm!=0){
    set.seed(seed)
    phenos[swarm=="A", id:=sample(id)]
    phenos[swarm=="B", id:=sample(id)]
  }
  phenos=phenos[match(rownames(loco.list[[1]]), phenos$id),]
  setnames(phenos, "id", "scanID")
  scanAnnot <- ScanAnnotationDataFrame(phenos[,.(scanID, temp.rack.cal, generation, photoperiod, diapause.bin, diapause.bin9, eggP, swarm)])
  
  filters=fread("/scratch/pae3g/final_reconstruction2/hwe_missing_maf_filters.txt")
  
  pass=filters[map_filter=="PASS", snp.id]
  rm(filters)

  nullmod.list=list()
  assoc.list=list()

  for (i in c("2L", "2R", "3L", "3R", "X")){

    chr.i=paste("chr", i, sep="")
    print(chr.i)
    if(pop!="both"){
    nullmod <- fitNullMM(scanData = scanAnnot, 
                       outcome = phenotype,  
                       covMatList = loco.list[[chr.i]], 
                       covars = c("temp.rack.cal", "generation", "photoperiod"),
                       family = binomial,
                       dropZeros=F)
    } else {
      nullmod <- fitNullMM(scanData = scanAnnot, 
                           outcome = phenotype,  
                           covMatList = loco.list[[chr.i]], 
                           covars = c("temp.rack.cal", "generation", "photoperiod", "swarm"),
                           family = binomial,
                           dropZeros=F)
      
    }
    nullmod.list[[i]]=nullmod
  
  ### do association w/ parallelization
    chrom <- getChromosome(geno)
    snpIDs <- getSnpID(geno, index=(chrom==i))
    
    snpIDs <- intersect(snpIDs, pass)
  
    groups <- split(snpIDs, ceiling(seq_along(snpIDs)/ceiling(length(snpIDs)/16)))
  
    assoc <- foreach(i=1:length(groups))%dopar%{
    
      model <- assocTestMM(genoData = genoData, 
                         nullMMobj = nullmod, 
                         test = "Score",
                         chromosome=i,
                         snp.include=groups[[i]])
      return(model)
    }
    assoc.list[[i]]=rbindlist(assoc)
  }

  assoc.results=rbindlist(assoc.list)
  assoc.results[,lod.geno:=-log10(Score.pval)]
  assoc.results[,pos:=getPosition(geno, snpID)]
  assoc.results[,perm:=perm]
  assoc.results[,seed:=seed]

  write.table(assoc.results, paste("/scratch/pae3g/final_reconstruction2/genesis_", phenotype, "~generation+temp.rack.cal+photoperiod_", pop, "_draw", draw, "_perm", perm, "_replaced.txt", sep=""), quote=FALSE, row.names=FALSE, sep="\t")


