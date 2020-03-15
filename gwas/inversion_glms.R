library(data.table)
library(cowplot)

#read in phenotypes and inversion calls based on Kapun 2013
phenos <- fread("/nv/vol186/bergland-lab/Priscilla/phenos_wolbachia.txt")
load("/scratch/pae3g/revisions/parents_hybrids_karyotype_calls.Rdata")

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

tableS3<-data.table(chromosome=c("3R", "3R", "3R", "2L", "2R"),
                    inversion=c("Mo", "C", "Payne", "t", "Ns"),
                   
                    
                    stage_8_both=c(summary(glm(diapause.bin~temp.rack.cal+photoperiod+swarm+generation+wolbachia+Mo, data=phenos, family="binomial"))$coefficients["Mo", "Pr(>|z|)"],
                                   summary(glm(diapause.bin~temp.rack.cal+photoperiod+swarm+generation+wolbachia+C, data=phenos, family="binomial"))$coefficients["C", "Pr(>|z|)"],
                                   summary(glm(diapause.bin~temp.rack.cal+photoperiod+swarm+generation+wolbachia+Payne, data=phenos[Payne!=2], family="binomial"))$coefficients["Payne", "Pr(>|z|)"],
                                   summary(glm(diapause.bin~temp.rack.cal+photoperiod+swarm+generation+wolbachia+t, data=phenos, family="binomial"))$coefficients["t", "Pr(>|z|)"],
                                   summary(glm(diapause.bin~temp.rack.cal+photoperiod+swarm+generation+wolbachia+Ns, data=phenos, family="binomial"))$coefficients["Ns", "Pr(>|z|)"]),
                    
                    stage_10_both=c(summary(glm(diapause.bin9~temp.rack.cal+photoperiod+swarm+generation+wolbachia+Mo, data=phenos, family="binomial"))$coefficients["Mo", "Pr(>|z|)"],
                                    summary(glm(diapause.bin9~temp.rack.cal+photoperiod+swarm+generation+wolbachia+C, data=phenos, family="binomial"))$coefficients["C", "Pr(>|z|)"],
                                    summary(glm(diapause.bin9~temp.rack.cal+photoperiod+swarm+generation+wolbachia+Payne, data=phenos[Payne!=2], family="binomial"))$coefficients["Payne", "Pr(>|z|)"],
                                    summary(glm(diapause.bin9~temp.rack.cal+photoperiod+swarm+generation+wolbachia+t, data=phenos, family="binomial"))$coefficients["t", "Pr(>|z|)"],
                                    summary(glm(diapause.bin9~temp.rack.cal+photoperiod+swarm+generation+wolbachia+Ns, data=phenos, family="binomial"))$coefficients["Ns", "Pr(>|z|)"]),
                    
                    
                    stage_8_A=c(summary(glm(diapause.bin~temp.rack.cal+photoperiod+generation+wolbachia+Mo, data=phenos[swarm=="A"], family="binomial"))$coefficients["Mo", "Pr(>|z|)"],
                                summary(glm(diapause.bin~temp.rack.cal+photoperiod+generation+wolbachia+C, data=phenos[swarm=="A"], family="binomial"))$coefficients["C", "Pr(>|z|)"],
                                summary(glm(diapause.bin~temp.rack.cal+photoperiod+generation+wolbachia+Payne, data=phenos[Payne!=2][swarm=="A"], family="binomial"))$coefficients["Payne", "Pr(>|z|)"],
                                summary(glm(diapause.bin~temp.rack.cal+photoperiod+generation+wolbachia+t, data=phenos[swarm=="A"], family="binomial"))$coefficients["t", "Pr(>|z|)"],
                                summary(glm(diapause.bin~temp.rack.cal+photoperiod+generation+wolbachia+Ns, data=phenos[swarm=="A"], family="binomial"))$coefficients["Ns", "Pr(>|z|)"]),
                    
                    
                    stage_10_A=c(summary(glm(diapause.bin9~temp.rack.cal+photoperiod+generation+wolbachia+Mo, data=phenos[swarm=="A"], family="binomial"))$coefficients["Mo", "Pr(>|z|)"],
                                 summary(glm(diapause.bin9~temp.rack.cal+photoperiod+generation+wolbachia+C, data=phenos[swarm=="A"], family="binomial"))$coefficients["C", "Pr(>|z|)"],
                                 summary(glm(diapause.bin9~temp.rack.cal+photoperiod+generation+wolbachia+Payne, data=phenos[Payne!=2][swarm=="A"], family="binomial"))$coefficients["Payne", "Pr(>|z|)"],
                                 summary(glm(diapause.bin9~temp.rack.cal+photoperiod+generation+wolbachia+t, data=phenos[swarm=="A"], family="binomial"))$coefficients["t", "Pr(>|z|)"],
                                 summary(glm(diapause.bin9~temp.rack.cal+photoperiod+generation+wolbachia+Ns, data=phenos[swarm=="A"], family="binomial"))$coefficients["Ns", "Pr(>|z|)"]),
                   
                    
                    stage_8_B=c(summary(glm(diapause.bin~temp.rack.cal+photoperiod+generation+wolbachia+Mo, data=phenos[swarm=="B"], family="binomial"))$coefficients["Mo", "Pr(>|z|)"],
                                NA,
                                summary(glm(diapause.bin~temp.rack.cal+photoperiod+generation+wolbachia+Payne, data=phenos[Payne!=2][swarm=="B"], family="binomial"))$coefficients["Payne", "Pr(>|z|)"],
                                summary(glm(diapause.bin~temp.rack.cal+photoperiod+generation+wolbachia+t, data=phenos[swarm=="B"], family="binomial"))$coefficients["t", "Pr(>|z|)"],
                                NA),
                    
                    stage_10_B=c(summary(glm(diapause.bin9~temp.rack.cal+photoperiod+generation+wolbachia+Mo, data=phenos[swarm=="B"], family="binomial"))$coefficients["Mo", "Pr(>|z|)"],
                                 NA,
                                 summary(glm(diapause.bin9~temp.rack.cal+photoperiod+generation+wolbachia+Payne, data=phenos[Payne!=2][swarm=="B"], family="binomial"))$coefficients["Payne", "Pr(>|z|)"],
                                 summary(glm(diapause.bin9~temp.rack.cal+photoperiod+generation+wolbachia+t, data=phenos[swarm=="B"], family="binomial"))$coefficients["t", "Pr(>|z|)"],
                                NA))
                    
write.csv(tableS3, "/scratch/pae3g/revisions/figures/tables3.csv")