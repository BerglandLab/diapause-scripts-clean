#install.packages("car")

library(car)
library(data.table)
p<-fread("~/Desktop/phenos_wolbachia.txt")
a<-glm(diapause.bin~temp.rack.cal*photoperiod+generation+swarm+wolbachia, data=p, family="binomial")

b<-glm(diapause.bin9~temp.rack.cal*photoperiod+generation+photoperiod+swarm+wolbachia, data=p, family="binomial")




y<-as.data.table(Anova(a,"II", test="F"))
setnames(y, "Sum Sq", "sumsq")
y[,pve:=sumsq/sum(sumsq)*100]

x<-as.data.table(Anova(b,"II", test="F"))
setnames(x, "Sum Sq", "sumsq")
x[,pve:=sumsq/sum(sumsq)*100]


results<-cbind(y, x)
write.csv(results, "~/Box Sync/scripts-working/revisions/phenotype_variance.csv")


fisher.test(table(p$wolbachia,p$diapause.bin9))
fisher.test(table(p$wolbachia,p$diapause.bin))

