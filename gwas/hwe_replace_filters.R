library(GWASTools)
library(data.table)
library(SNPRelate)
library(gdsfmt)
library(doMC)
registerDoMC(20)
library(foreach)

##this script generates a master data table with informtion about each snp that is used for filtering and for permutation analysis, including:
    # of times the snp was heterozygous in a founder
    # of times the snp was replaced due to being missing in parent or in a "bad path" (haplotype path < 1 Mb)
    #hardy weinberg equilibrium of each snp in cage A, cage B, and combined
    #fst between cage A and cage B
    # population allele frequency

#calculate the number of heterozygous parental sites in each individual
bad.paths=fread("/scratch/pae3g/final_reconstruction2/bad_paths.csv")


import_data_het<-function(fn){
    sample=strsplit(fn,split="[.]")[[1]][1]
    if (!file.size(fn) == 0) {
        #read in data and rename columns
        dat<-fread(paste('zcat ', fn, sep=""))
        names(dat)=c("chr", "position", "h1", "h2")
        dat[,sample.id:=sample]
        print(sample)
        setkey(dat, chr, position)
        bad.paths.dat=bad.paths[sample.id==sample]
        #for each bad path segment, mask the appropriate haplotype
        foreach(bad.chr=bad.paths.dat$chr, bad.start=bad.paths.dat$cons.start, bad.stop=bad.paths.dat$cons.stop, bad.hap=bad.paths.dat$haplotype)%do%{
            if(bad.hap=="par1"){
                dat[chr==bad.chr & position>=bad.start & position<=bad.stop, h1:="."]
            }
            if(bad.hap=="par2") {
                dat[chr==bad.chr & position>=bad.start & position<=bad.stop, h2:="."]
            }
        }
        return(rbind(dat[h1=="0/1", .(chr, position, sample.id)], dat[h2=="0/1", .(chr, position, sample.id)]))
    }
}

swarm.het<-foreach(swarm=c("A", "B"))%do%{
    
    wd=paste("/scratch/pae3g/genome-reconstruction/", swarm,  sep="")
    setwd(wd)
    file_list <- list.files(pattern = "PAE.*.14.vcf.gz$")
    obj <- foreach(fn=file_list)%dopar%{
        import_data_het(fn)
    }
    obj=rbindlist(obj)
    return(obj)
}

swarm.het=rbindlist(swarm.het)

swarm.het.sum=swarm.het[,.(n=.N), .(chr, position)][order(chr, position)]
swarm.het.samp=swarm.het[,.(n=.N), .(chr, position, sample.id)][order(chr, position)]

write.table(swarm.het.sum, "/scratch/pae3g/final_reconstruction2/parent_het_per_snp.txt", sep="\t", quote=F, row.names=F)


#calculate HWE for every SNP in swarm A, swarm B, and both 
geno <- snpgdsOpen("/scratch/pae3g/genome-reconstruction/final2_draw1_replaced.vcf.gds", allow.fork=TRUE, readonly = FALSE)
phenos <- fread("/scratch/pae3g/phenos/phenos_062018.txt")
setkey(phenos, id)
genotyped_samples<- read.gdsn(index.gdsn(geno, "sample.id"))


#add in metadata to gds to file
metadata=fread("/scratch/pae3g/hybrid/PAE_AOB_library_metadata_all.txt", header=FALSE)
names(metadata)=c("sample.id", "generation", "library", "swarm")
metadata=metadata[sample.id%in%genotyped_samples]
metadata=metadata[match(genotyped_samples, metadata$sample.id),]

samp.annot <- metadata[,.(generation, swarm)]
add.gdsn(geno, "sample.annot", samp.annot)


#make list of A and B samples
A=intersect(phenos[swarm=="A", id],genotyped_samples)
B=intersect(phenos[swarm=="B", id],genotyped_samples)

#calculate HWE for each swarm and combine with snp information
a.hwe=snpgdsHWE(geno, sample.id=A, snp.id=NULL, with.id=T)
a.info=snpgdsSNPList(geno, sample.id=A)
a.dat=data.table(snp.id=a.hwe$snp.id,id=a.info$rs.id,chr=a.info$chromosome, pos=a.info$position, a.hwe=a.hwe$pvalue, a.freq=a.info$afreq)

b.hwe=snpgdsHWE(geno, sample.id=B, snp.id=NULL, with.id=T)
b.info=snpgdsSNPList(geno, sample.id=B)
b.dat=data.table(snp.id=b.hwe$snp.id, id=b.info$rs.id,chr=b.info$chromosome, pos=b.info$position, b.hwe=b.hwe$pvalue, b.freq=b.info$afreq)

hwe=snpgdsHWE(geno,  snp.id=NULL, with.id=T)
info=snpgdsSNPList(geno )
dat=data.table(snp.id=hwe$snp.id, id=info$rs.id,chr=info$chromosome, pos=info$position, hwe=hwe$pvalue, freq=info$afreq)

#merge all together
hwe.filter=merge(a.dat, merge(dat, b.dat, by=c("snp.id", "id", "chr", "pos")), by=c("snp.id", "id", "chr", "pos"))

#load in genotype replacement data from the logs produced by replace_missing_data_vcf.sh

a.replace=fread("/scratch/pae3g/genome-reconstruction/replace_A_1_better.log")
a.replace=a.replace[2:.N] #removes #CHROM row
names(a.replace)=c("chr", "pos", "nRef.a", "nAlt.a", "mostLikelyGeno.a", "n.a", "total.a")
b.replace=fread("/scratch/pae3g/genome-reconstruction/replace_B_1_better.log")
b.replace=b.replace[2:.N, 3:7] #removes #CHROM row
names(b.replace)=c("nRef.b", "nAlt.b", "mostLikelyGeno.b", "n.b", "total.b")

#merge with other data
replace=cbind(a.replace, b.replace)
replace[,pos:=as.integer(pos)]

hwe.filter=merge(hwe.filter, replace, by=c("chr", "pos"))
hwe.filter=hwe.filter[chr%in%(c("2L", "2R", "3L", "3R", "X"))]


# Fst estimation between cages A and B

group <- as.factor(read.gdsn(index.gdsn(geno, "sample.annot/swarm")))
v <- snpgdsFst(geno, population=group, method="W&C84",autosome.only=FALSE, with.id=T, remove.monosnp=F)

fst=data.table(snp.id=v$snp.id, fst=v$FstSNP)

#merge fst with other data
hwe.filter=merge(hwe.filter, fst, by="snp.id")

#merge in heterozygosity data calculated aboave
het=fread("/scratch/pae3g/final_reconstruction2/parent_het_per_snp.txt")
setnames(het, "position", "pos")
setnames(het, "n", "n.het")

hwe.filter=merge(hwe.filter, het, by=c("chr", "pos"), all=T)
hwe.filter[is.na(n.het), n.het:=0]
#calculate the total number of individuals that had imputed (heterozygous parents)
hwe.filter[,n.imputed:=n.a+n.b]

#calculate the total number of repacedsites 
hwe.filter[,total.replaced:=n.a+n.b+n.het]

#rename as "dat"
dat<-hwe.filter


#make a filter column
dat[,qc_filter:="PASS"]
dat[,map_filter:="PASS"]


#filter snps with more than 10% of genotypes replaced

dat[total.replaced>566, qc_filter:="FAIL"]

#filter things out of hwe in one or both cages
dat[a.hwe<10^-20, qc_filter:="FAIL"] #12800
dat[b.hwe<10^-20, qc_filter:="FAIL"] #13000
dat[hwe<10^-20, qc_filter:="FAIL"] #16000

#filter things with high fst between pops
dat[fst>.2, qc_filter:="FAIL"] #13000

#filter things that are monomorphic within a cage
dat[a.freq==0|a.freq==1|b.freq==0|b.freq==1, qc_filter:="FAIL"] #~1 million sites

#leaves ~1.2 million snps 

#for mapping don't have to be as stringent; can keep things that are monomorphic in a cage
dat[total.replaced>566, map_filter:="FAIL"]

#filter things out of hwe in one or both cages
dat[a.hwe<10^-20, map_filter:="FAIL"] #12800
dat[b.hwe<10^-20, map_filter:="FAIL"] #13000
dat[hwe<10^-20, map_filter:="FAIL"] #16000
dat[fst>.2, map_filter:="FAIL"] 
dat[freq==0, map_filter:="FAIL"]

#this should result in a good snpset for PCA and gwas to save computational time

write.table(dat, "/scratch/pae3g/final_reconstruction2/hwe_missing_maf_filters.txt", sep="\t", quote=F, row.names=F)


#add in bin information and create new filter for adaptive permutations

library(data.table)
library(Hmisc)

dat<-fread("/scratch/pae3g/final_reconstruction2/hwe_missing_maf_filters.txt")
dat[,snpID:=snp.id]

dat[,MAF:=ifelse(freq>.5, 1-freq, freq)]
dat[,mis.bin:=cut2(total.replaced, cuts=seq(0,600, by=100))]
dat[,maf.bin:=cut2(MAF, cuts=seq(0, .5, by=.05))]

bins<-dat[map_filter=="PASS", .(n.in.bin=.N), .(maf.bin, mis.bin)]

dat<-merge(dat, bins, all=T, by=c("maf.bin", "mis.bin"))

dat[, adaptive_perm_filter:=map_filter]
dat[n.in.bin>10000, adaptive_perm_filter:="FAIL"]

write.table(dat, "/scratch/pae3g/final_reconstruction2/hwe_missing_maf_filters.txt", sep="\t", quote=F, row.names=F)

