#get list of genotyped sample names:
library(data.table)
library(foreach)


#list the files
tim.geno<-foreach(swarm=c("A", "B"), .combine="rbind", .errorhandling="remove")%do%{
    setwd(paste0("/scratch/pae3g/genome-reconstruction/", swarm))
    fn=list.files(pattern="PAE.*estimate.14.haps")
    dat<-foreach(f=fn, .combine="rbind", .errorhandling="remove")%do%{
        h<-fread(f)
        #print(h)
        samp=tstrsplit(f, split="[.]")[[1]]
        print(samp)
        p1<-h[chromosome=="2L"&(stop>=3504474)&(start<=3504474), par1]
        #print(p1)
        p2<-h[chromosome=="2L"&(stop>=3504474)&(start<=3504474), par2]
        #print(p2)
        return(data.table(id=samp,
                          par1=p1,
                          par2=p2,
                          swarm=swarm))
    }
    return(dat)
}

key<-fread("/scratch/pae3g/hybrid/parent_line_by_swarm.txt", header=F)
setnames(key, c("sample.id", "num.68", "num.34", "swarm"))

key[,par1:=num.34]
key[,par2:=num.34]
key[,line1:=sample.id]
key[,line2:=sample.id]

tim.geno<-merge(tim.geno, key[,.(par1, line1, swarm)], by=c("par1", "swarm"))
tim.geno<-merge(tim.geno, key[,.(par2, line2, swarm)], by=c("par2", "swarm"))

#save this in leased storage and switch to workstation

write.csv(tim.geno, "/nv/vol186/bergland-lab/Priscilla/tim_haplotypes.csv", row.names=F)


#workstation
library(SNPRelate)
library(data.table)
library(cowplot)
snpgdsVCF2GDS(vcf.fn = "/mnt/internal_2/priscilla/HSparents/variants/merged.hc.raw.vcf",out.fn = "/mnt/internal_2/priscilla/HSparents/variants/merged.hc.raw.vcf.gds", method="copy.num.of.ref", verbose=T,)
#pull in full vcf

geno<-snpgdsOpen("/mnt/internal_2/priscilla/HSparents/variants/merged.hc.raw.vcf.gds")

#get variant id of lstim snp (listed at 2L:3504474)

a<-snpgdsSNPList(gdsobj=geno)
samps<-data.table(snp.id=a$snp.id,
                  chr=a$chromosome,
                  pos=a$pos)

samps[pos==3504474]
y<-snpgdsGetGeno(geno, snp.id=167773, with.id=T)

haps<-data.table(sample.id=y$sample.id,
                 geno1=y$genotype[,1],
                 geno2=y$genotype[,1],
                 line1=y$sample.id,
                 line2=y$sample.id)

write.csv(haps, "/mnt/sammas_storage/bergland-lab/Priscilla/tim_parental_genos.csv", row.names=F)


tim.geno<-merge(tim.geno, haps[,.(geno1,line1)],  "line1")
tim.geno<-merge(tim.geno, haps[,.(geno2, line2)], "line2")

tim.geno[geno1==0&geno2==0, final.geno:=0]
tim.geno[geno1==0&geno2==2, final.geno:=1]
tim.geno[geno1==2&geno2==0, final.geno:=1]
tim.geno[geno1==2&geno2==2, final.geno:=2]

write.csv(tim.geno, "/scratch/pae3g/tim_genotypes.csv")

