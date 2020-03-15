library(data.table)
library(SNPRelate)
library(doMC)
registerDoMC(20)

#i=1

phenos <- fread("/scratch/pae3g/oldscratch_recovered/phenos/phenos_062018.txt")
filters=fread("/scratch/pae3g/oldscratch_recovered/final_reconstruction2/hwe_missing_maf_filters.txt")
pass=filters[map_filter=="PASS", snp.id]

foreach(i=c(4:100))%do%{
    filestem<- paste0("/scratch/pae3g/genome-reconstruction/final2_draw", i, "_replaced")
    genofile<-snpgdsOpen(paste0(filestem,".vcf.gds"))
    genotyped_samples<- read.gdsn(index.gdsn(genofile, "sample.id"))
    
    ids.to.use=intersect(genotyped_samples, phenos$sample.id)
    
    pca.obj <- snpgdsPCA(gdsobj=genofile,
                         sample.id=genotyped_samples,
                         snp.id=pass,
                         autosome.only=F,
                         maf=0.01,
                         num.thread=20)

    save(pca.obj, file=paste0(filestem, "_PCA_maf0.01.Rdata"))
}