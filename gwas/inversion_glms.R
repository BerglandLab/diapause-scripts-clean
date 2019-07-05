library(data.table)
library(cowplot)

#read in phenotypes and inversion calls based on Kapun 2013
phenos <- fread("/scratch/pae3g/phenos/phenos_062018.txt")
kary.ag<-fread("/scratch/pae3g/evolution/final3.vcf.karytypecalls.csv", header=T)

phenos<-merge(kary.ag, phenos, by="sample.id")


#pull out individual inversions in 0/1/2 format

phenos[,C:=0]
phenos[,Mo:=0]
phenos[,Payne:=0]
phenos[chr3R=="In(3R)Payne;In(3R)Mo" , Mo:=1]
phenos[chr3R=="In(3R)Payne;In(3R)Mo" , Payne:=1]
phenos[chr3R=="std;In(3R)Mo", Mo:=1]
phenos[chr3R=="std;In(3R)C", C:=1]
phenos[chr3R=="std;In(3R)Payne", Payne:=1]
phenos[chr3R=="In(3R)Payne;In(3R)C", Payne:=1]
phenos[chr3R=="In(3R)Payne;In(3R)C", C:=1]
phenos[chr3R=="In(3R)C", C:=2]
phenos[chr3R=="In(3R)Payne", Payne:=2]
phenos[chr3R=="In(3R)Mo", Mo:=2]
phenos[, Ns:=0]
phenos[chr2R=="In(2R)Ns", Ns:=2]
phenos[chr2R=="chr2Rstd;In(2R)Ns", Ns:=1]
phenos[, t:=0]
phenos[chr2L=="In(2L)t", t:=2]
phenos[chr2L=="chr2Lstd;In(2L)t", t:=1]

#test individual inversions
summary(glm(diapause.bin9~temp.rack.cal+photoperiod+swarm+generation+Mo, data=phenos, family="binomial"))
summary(glm(diapause.bin9~temp.rack.cal+photoperiod+swarm+generation+C, data=phenos, family="binomial"))
summary(glm(diapause.bin9~temp.rack.cal+photoperiod+swarm+generation+Payne, data=phenos[Payne!=2], family="binomial"))
summary(glm(diapause.bin9~temp.rack.cal+photoperiod+swarm+generation+t, data=phenos, family="binomial"))
summary(glm(diapause.bin9~temp.rack.cal+photoperiod+swarm+generation+Ns, data=phenos, family="binomial"))
#all chr3R genos together
summary(glm(diapause.bin9~temp.rack.cal+photoperiod+swarm+generation+chr3R, data=phenos, family="binomial"))


