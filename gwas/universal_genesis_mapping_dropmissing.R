#!/usr/bin/env Rscript

library(GWASTools, quietly=TRUE)
library(GENESIS, quietly=TRUE)
library(data.table, quietly=TRUE)
library(SNPRelate, quietly=TRUE)
library(gdsfmt, quietly=TRUE)
library(doMC, quietly=TRUE)
registerDoMC(20)
library(foreach, quietly=TRUE)
library(biglasso)


args=commandArgs(trailingOnly=TRUE)

pop=args[1]
phenotype=args[2]
draw=args[3]
perm=args[4]
seed=args[5]
model=args[6]


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
#i=1

filestem<- paste0("/scratch/pae3g/genome-reconstruction/final2_draw", draw, "_replaced")

#only run if file doesn't exist



### genoData 	
geno <- GdsGenotypeReader.par(filename = paste(filestem, ".vcf.gds", sep=""), allow.fork=T)

g <- GenotypeData(geno)

#read data
phenos <- fread("/scratch/pae3g/revisions/phenos_wolbachia_missinggeno.txt")
setkey(phenos, id)

### load appropriate GRM(s)	
phenos[,sex := "F"]	
if (pop=="both"&model=="nonloco"){
  load(paste(filestem, "_nonLOCO_Eigenstrat_allgwasSNPs.Rdat", sep=""))
}else if (pop!="both"&model=="nonloco") {
  load(paste(filestem, "_nonLOCO_Eigenstrat_allgwasSNPs_", pop, ".Rdat", sep=""))
}else if (pop!="both"&model=="loco"){
  load(paste(filestem, "_LOCO_Eigenstrat_allsnps_", pop, ".Rdat", sep=""))
  a<-loco.list[[1]]
}else if (pop=="both"&model=="loco"){
  load(paste(filestem, "_LOCO_Eigenstrat_allsnps.Rdat", sep=""))
  a<-loco.list[[1]]
}

#drop samples with reconstruction issues

phenos <- phenos[n.chr.imp==5&prop.unknown<=0.05]

#for single swarm analysis, subset phenos to only that swarm
if(pop!="both"){
  phenos<-phenos[swarm==pop]
}

#subset sample GRM to only the samples that will be used
a<-a[phenos$sample.id, phenos$sample.id]

#if a permutation, randomly scramble sample ids within a swarm. if perm is 0, ids will stay the same
if(perm!=0){
  set.seed(seed)
  phenos[swarm=="A", id:=sample(id)]
  phenos[swarm=="B", id:=sample(id)]
}
phenos=phenos[match(rownames(a), phenos$id),]
phenos[,scanID:=id]

#set mapping phenotype
if(phenotype=="diapause.bin"){
  phenos[,mapping_pheno:=diapause.bin]
}
if(phenotype=="diapause.bin9"){
  phenos[,mapping_pheno:=diapause.bin9]
}
if(!file.exists(paste("/scratch/pae3g/revisions/genesis_", phenotype, "_draw", draw, "_perm", perm, "_pop", pop, "_" , model, "_allsnpgrm_wolbachia_dropmissing.Rdat", sep=""))){  
  print("starting gwas")
  #make scanAnnot dataframe
  scanAnnot <- ScanAnnotationDataFrame(phenos[,.(scanID, temp.rack.cal, generation, photoperiod, wolbachia, mapping_pheno, swarm)])
  
  #filter down to only mapping snps
  filters=fread("/scratch/pae3g/oldscratch_recovered/final_reconstruction2/hwe_missing_maf_filters.txt")
  
  pass=filters[map_filter=="PASS", snp.id]
  #rm(filters)
  
  #nonloco mapping
  
  if(model=="nonloco"){
    
    #include swarm as covariate if both
    if(pop=="both"){
      nullmod <- fitNullModel(x = scanAnnot, 
                              outcome = "mapping_pheno", 
                              cov.mat = a, 
                              covars = c("temp.rack.cal", "generation", "photoperiod", "swarm", "wolbachia"),
                              family = binomial,
                              drop.zeros=F)
    }else {
      nullmod <- fitNullModel(x = scanAnnot, 
                              outcome = "mapping_pheno",  
                              cov.mat = a, 
                              covars = c("temp.rack.cal", "generation", "photoperiod", "wolbachia"),
                              family = binomial,
                              drop.zeros=F)
    }
    
    
    ### do association w/ parallelization
    snpIDs <- getSnpID(geno)
    
    snpIDs <- intersect(snpIDs, pass)
    
    groups <- split(snpIDs, ceiling(seq_along(snpIDs)/ceiling(length(snpIDs)/20)))
    assoc <- foreach(i=1:length(groups))%dopar%{
      genoIterator <- GenotypeBlockIterator(g, snpInclude=groups[[i]])
      
      genesis.model <- assocTestSingle(genoIterator, 
                                       null.model = nullmod, 
                                       test = "Score",
                                       verbose=T)
      # snp.include=groups[[i]])
      return(genesis.model)
    }
    
    
    
    assoc.results=rbindlist(assoc)
    
  }
  
  
  #loco mapping
  if(model=="loco"){
    assoc.list=list()
    for (i in c("2L", "2R", "3L", "3R", "X")){
      
      chr.i=paste("chr", i, sep="")
      pass=filters[map_filter=="PASS"&chr==i, snp.id]
      
      print(chr.i)
      if(pop!="both"){
        #make a new null model at the beginning of each chromosome
        if(i=="2L"|i=="3L"|i=="X"){
          nullmod <- fitNullModel(x = scanAnnot, 
                                  outcome = "mapping_pheno",  
                                  cov.mat = loco.list[[chr.i]][phenos$sample.id, phenos$sample.id], 
                                  covars = c("temp.rack.cal", "generation", "photoperiod", "wolbachia"),
                                  family = binomial,
                                  drop.zeros=F)
        }
      } else {
        if(i=="2L"|i=="3L"|i=="X"){
          nullmod <- fitNullModel(x = scanAnnot, 
                                  outcome = "mapping_pheno",  
                                  cov.mat = loco.list[[chr.i]][phenos$sample.id, phenos$sample.id], 
                                  covars = c("temp.rack.cal", "generation", "photoperiod", "swarm", "wolbachia"),
                                  family = binomial,
                                  drop.zeros=F)
        }
      }
      
      
      ### do association w/ parallelization
      snpIDs <- getSnpID(geno)
      
      snpIDs <- intersect(snpIDs, pass)
      
      groups <- split(snpIDs, ceiling(seq_along(snpIDs)/ceiling(length(snpIDs)/20)))
      assoc <- foreach(j=1:length(groups))%dopar%{
        genoIterator <- GenotypeBlockIterator(g, snpInclude=groups[[j]])
        
        genesis.model <- assocTestSingle(genoIterator, 
                                         null.model = nullmod, 
                                         test = "Score",
                                         verbose=T)
        # snp.include=groups[[i]])
        
        return(genesis.model)
      }
      
      assoc.list[[i]]=rbindlist(assoc)
    }
    
    
    assoc.results=rbindlist(assoc.list)
    
  }
  assoc.results[,lod.geno:=-log10(Score.pval)]
  assoc.results[,pos:=getPosition(geno, variant.id)]
  assoc.results[,chr:=getChromosome(geno, variant.id)]
  
  assoc.results[,perm:=perm]
  assoc.results[,seed:=seed]
  assoc.results[,model:=model]
  assoc.results[,phenotype:=phenotype]
  
  save(assoc.results, file=paste("/scratch/pae3g/revisions/genesis_", phenotype, "_draw", draw, "_perm", perm, "_pop", pop, "_" , model, "_allsnpgrm_wolbachia_dropmissing.Rdat", sep=""))
}else{
  print("gwas already done; loading results")
  load(paste("/scratch/pae3g/revisions/genesis_", phenotype, "_draw", draw, "_perm", perm, "_pop", pop, "_" , model, "_allsnpgrm_wolbachia_dropmissing.Rdat", sep=""))
}



##### run LASSO #######
outdir<-"/scratch/pae3g/revisions/lasso/"

if(!file.exists(paste0(outdir, "lasso_", phenotype, "_draw",draw, "_perm", perm, "_", model, "_pop", pop, ".dropmissing.lassoSites.Rdata"))){
  print("running lasso")
  maf.threshold=0.05
  maxRank <- 10000   ### how many SNPs to dump into LASSO?
  
  #basic info about snps
  info<-fread("/scratch/pae3g/revisions/snpinfo.txt")
  
  #phenotypes --treat identically to how they were treated in GENESIS so that the permutation happens the same way
  setkey(phenos, id)
  
  
  pheno.small <- phenos[,c("id", "mapping_pheno", "temp.rack.cal", "swarm", "generation", "photoperiod", "wolbachia"), with=F]
  setnames(pheno.small, c("temp.rack.cal", "generation", "photoperiod"), c("temp", "gen", "photo"))
  pheno.small[,orig_id:=id]
  pheno.small <- pheno.small[!is.na(mapping_pheno)]
  
  
  #load PCA
  print("loading PCA")
  load(paste0(filestem, "_PCA_maf0.01.Rdata"))
  
  pca.dt <- as.data.table(pca.obj$eigenvect)
  setnames(pca.dt, names(pca.dt), gsub("V", "PC_", names(pca.dt)))
  pca.dt[,id:=pca.obj$sample.id]
  
  
  gwas <- assoc.results
  
  #filter gwas for MAF
  gwas<-gwas[freq>=maf.threshold&freq<(1-maf.threshold)]
  gwas[,rs:=paste(chr, pos, "SNP", sep="_")]
  
  #rank GWAS
  gwas[,r:=frank(Score.pval, na.last=T, ties.method="min")]
  
  #make genotype matrix for top maxRank SNPs
  
  close(geno)
  geno<-snpgdsOpen( paste(filestem, ".vcf.gds", sep=""))
  print('making genomat')
  genomat<-as.data.table(snpgdsGetGeno(geno, snp.id=gwas[r<maxRank, variant.id], with.id=FALSE, verbose=TRUE))
  
  setnames(genomat,
           names(genomat),
           gwas[r<maxRank, rs])
  genomat[,id:=read.gdsn(index.gdsn(geno, "sample.id"))]
  setkey(genomat, id)
  
  Y <- pheno.small[id%in%genomat$id, c("mapping_pheno", "id"), with=F]
  setkey(Y, id)
  
  print("generating models")
  X.bm <- list()
  
  ## environment only
  ### extract out
  setkey(pheno.small, id)
  pheno.small.tmp <- pheno.small[J(Y$id)]
  pheno.small.tmp[,swarm:=as.numeric(as.factor(swarm))]
  pheno.small.tmp[,wolbachia:=as.numeric(as.factor(wolbachia))]
  
  
  ### center
  pheno.small.tmp[,temp.c := temp - mean(temp)]
  pheno.small.tmp[,gen.c := gen - mean(gen)]
  pheno.small.tmp[,photo.c := gen - mean(photo)]
  pheno.small.tmp[,swarm.c := gen - mean(swarm)]
  pheno.small.tmp[,wolbachia.c := gen - mean(wolbachia)]
  
  
  ### make big.matrix
  X.bm$env <- as.big.matrix(pheno.small.tmp[, c("temp.c", "swarm.c", "gen.c", "photo.c", "wolbachia.c"), with=F])
  
  ## environment + PC
  ### merge pheno + pc
  pheno.small.pc.tmp <- merge(pheno.small.tmp, pca.dt, by="id")
  
  #table(pheno.small.pc.tmp$id==pheno.small.tmp$id)
  
  ### make big.matrix
  X.bm$env.pc <- as.big.matrix(pheno.small.pc.tmp[, c("temp.c", "swarm.c", "gen.c", "photo.c", "wolbachia.c",
                                                      names(pheno.small.pc.tmp)[grepl("PC", names(pheno.small.pc.tmp))]), with=F])
  
  ## environment + PC + genotype
  ### merge pheno + pc + genotype
  setkey(genomat, id)
  pheno.small.pc.geno.tmp <- merge(pheno.small.pc.tmp, genomat)
  
  #table(pheno.small.pc.geno.tmp$id==pheno.small.tmp$id)
  
  ### make big.matrix
  X.bm$env.pc.geno <- as.big.matrix(pheno.small.pc.geno.tmp[, c("temp.c", "swarm.c", "gen.c", "photo.c", "wolbachia.c",
                                                                names(pheno.small.pc.geno.tmp)[grepl("PC", names(pheno.small.pc.geno.tmp))],
                                                                names(pheno.small.pc.geno.tmp)[grepl("SNP", names(pheno.small.pc.geno.tmp))]),
                                                            with=F])
  
  print("running lasso")
  
  cvfit.out <- foreach(i=1:length(X.bm))%do%{
    cv.biglasso(X=X.bm[[i]], y=Y$mapping_pheno, nfolds = 10, ncores = 10, family="binomial", screen="SSR-Slores", penalty="enet", max.iter=50000, alpha=.5)
  }
  
  
  sites <- foreach(i = 1:length(cvfit.out), .combine="rbind")%do%{
    # i <-2
    cvi <- cvfit.out[[i]]
    
    coefs.tmp <- coef(cvi)
    
    data.table(coef=coefs.tmp[which(coefs.tmp != 0)],
               site=dimnames(coefs.tmp)[[1]][which(coefs.tmp != 0)],
               lassoMod=names(X.bm)[i],
               maf=maf.threshold,
               perm=perm,
               loco=model,
               GRM=draw,
               pheno=phenotype,
               pop=pop)
    
  }
  
  pred <- foreach(i = 1:length(cvfit.out), .combine="rbind")%do%{
    cvi <- cvfit.out[[i]]
    
    
    data.table(pred=predict(cvi, X.bm[[i]])[,1],
               obs=Y$mapping_pheno,
               lassoMod=names(X.bm)[i],
               maf=maf.threshold,
               perm=perm,
               loco=model,
               GRM=draw,
               pheno=phenotype,
               pop=pop,
               id=Y$id)
    
    
  }
  
  print("saving output")
  pred.out.fn <- paste0(outdir, "lasso_", phenotype, "_draw",draw, "_perm", perm, "_", model, "_pop", pop, ".dropmissing.lassoPred.Rdata")
  
  sites.out.fn <-   paste0(outdir, "lasso_", phenotype, "_draw",draw, "_perm", perm, "_", model, "_pop", pop, ".dropmissing.lassoSites.Rdata")
  
  fits.out.fn <-  paste0(outdir, "lasso_", phenotype, "_draw",draw, "_perm", perm, "_", model, "_pop", pop, ".dropmissing.lassoData.Rdata")
  
  save(pred, file=pred.out.fn)
  
  save(sites, file=sites.out.fn)
  
  save(cvfit.out, file=fits.out.fn)
  
}

if(!file.exists(paste0(outdir, "GIF_", phenotype, "_draw",draw, "_perm", perm, "_", model, "_pop", pop, ".dropmissing.Rdata"))){
  #CALCULATE GIF and SAVE
  gwas<-assoc.results
  print("calculating GIF")
  chisq.05 <- qchisq(1-gwas[freq>=0.05 & freq <=0.95, Score.pval],1)
  lambda.05 = median(chisq.05)/qchisq(0.5,1)
  chisq.01<-qchisq(1-gwas[freq>=0.01& freq <=0.99, Score.pval],1)
  lambda.01 = median(chisq.01)/qchisq(0.5,1)
  g<-data.table(perm=perm,
                draw=draw,
                lambda.05=lambda.05,
                lambda.01=lambda.01,
                pop=pop,
                pheno=phenotype,
                model=model)
  
  save(g, file=paste0(outdir, "GIF_", phenotype, "_draw",draw, "_perm", perm, "_", model, "_pop", pop, ".dropmissing.Rdata"))
}

if(!file.exists(paste0(outdir, "GIFbychr_", phenotype, "_draw",draw, "_perm", perm, "_", model, "_pop", pop, ".dropmissing.Rdata"))){
  #CALCULATE GIF and SAVE
  x<-foreach(CHR=c("2L", "2R", "3L", "3R", "X"))%do%{
    gwas<-assoc.results[chr==CHR]
    print("calculating GIF")
    chisq.05 <- qchisq(1-gwas[freq>=0.05 & freq <=0.95, Score.pval],1)
    lambda.05 = median(chisq.05)/qchisq(0.5,1)
    chisq.01<-qchisq(1-gwas[freq>=0.01& freq <=0.99, Score.pval],1)
    lambda.01 = median(chisq.01)/qchisq(0.5,1)
    g<-data.table(perm=perm,
                  draw=draw,
                  lambda.05=lambda.05,
                  lambda.01=lambda.01,
                  pop=pop,
                  pheno=phenotype,
                  model=model,
                  chr=CHR)
    return(g)
  }
  x<-rbindlist(x)
  
  save(x, file=paste0(outdir, "GIFbychr_", phenotype, "_draw",draw, "_perm", perm, "_", model, "_pop", pop, ".dropmissing.Rdata"))
}
