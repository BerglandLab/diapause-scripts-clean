#!/usr/bin/env Rscript

library(GWASTools, quietly=TRUE)
library(GENESIS, quietly=TRUE)
library(data.table, quietly=TRUE)
library(SNPRelate, quietly=TRUE)
library(gdsfmt, quietly=TRUE)
library(doMC, quietly=TRUE)
registerDoMC(20)
library(foreach, quietly=TRUE)

args=commandArgs(trailingOnly=TRUE)

pop=args[1]
phenotype=args[2]
perm=args[3]
seed=args[4]

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

  #read data

filters=fread("/scratch/pae3g/final_reconstruction2/hwe_missing_maf_filters.txt")


draws<-foreach(draw=c(1:100))%do%{
  phenos <- fread("/scratch/pae3g/phenos/phenos_062018.txt")
  setkey(phenos, id)
  ### some twiddles	
  phenos[,sex := "F"]	
  if (pop=="both"){
    load(paste("/scratch/pae3g/genome-reconstruction/final2_draw", draw, "_replaced_LOCO_Eigenstrat_strict_both.Rdat", sep=""))
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
  
  
  pass=filters[adaptive_perm_filter=="PASS", snp.id]

  geno <- GdsGenotypeReader.par(filename = paste("/scratch/pae3g/genome-reconstruction/final2_draw", draw, "_replaced.vcf.gds", sep=""), allow.fork=T)
  genoData <- GenotypeData(geno)
  
  assoc.list<-foreach(i=c("2", "3", "X"))%do%{
    if(i=="2"){
      chr.i="chr2L" }   
    if(i=="3"){
      chr.i="chr3L"}
    if(i=="X"){
      chr.i="chrX"}
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

  ### do association w/ parallelization
    chrom <- getChromosome(geno)
    if(i=="2"){
      snpIDs <- getSnpID(geno, index=(chrom=="2L"|chrom=="2R"))}
    if(i=="3"){
      snpIDs <- getSnpID(geno, index=(chrom=="3L"|chrom=="3R"))}
    if(i=="X"){
      snpIDs <- getSnpID(geno, index=(chrom=="X"))}
    
    
    snpIDs <- intersect(snpIDs, pass)
  
    groups <- split(snpIDs, ceiling(seq_along(snpIDs)/ceiling(length(snpIDs)/20)))
  
    assoc <- foreach(i=1:length(groups))%dopar%{
    
      model <- assocTestMM(genoData = genoData, 
                         nullMMobj = nullmod, 
                         test = "Score",
                         snp.include=groups[[i]])
      return(as.data.table(model))
    }
    return(rbindlist(assoc))
  }

  assoc.results=rbindlist(assoc.list)
  assoc.results[,lod.geno:=-log10(Score.pval)]
  assoc.results[,pos:=getPosition(geno, snpID)]
  assoc.results[,perm:=perm]
  assoc.results[,seed:=seed]
  assoc.results[,draw:=draw]
  close(geno)
  return(assoc.results)
}
draws<-rbindlist(draws)

draws.sum=draws[,.(avg.pval=mean(Score.pval,na.rm=T), med.pval=median(Score.pval,na.rm=T), min.pval=min(Score.pval, na.rm=T), max.pval=max(Score.pval, na.rm=T), sd.pval=sd(Score.pval, na.rm=T)), .(snpID, chr, pos)]

draws.sum[,perm:=perm]
draws.sum[,seed:=seed]


 write.table(draws.sum, paste("/scratch/pae3g/final_reconstruction2/adaptive_perm_", perm, "_", pheno, ".txt", sep=""), quote=FALSE, row.names=FALSE, sep="\t")


