library(data.table)
library(foreach)
library(cowplot)
library(doMC)
registerDoMC(10)
load(file="/mnt/spicy_1/pigmentation/inputData/dat.Rdata")
setkey(dat, chr, pos, pop)
ids=as.data.table(read.csv("/mnt/pricey_2/priscilla/hybrid/all_popinfo_formap_set.csv"))
setnames(ids, "name", "pop")
setkey(ids, pop)

dat<-merge(dat, ids, by="pop")

#subset to core20
core=dat[set=="Core20"]

save(core, file="/mnt/pricey_2/priscilla/core20.rdat")
#calculate ref freq
core[,rf:=1-af]
#logit transform ref freq
core[,logit.rf:=qlogis(rf)]
core[,pop.char:=as.character(pop_name)]

setkey(core, pop.char)
deltas<-foreach(p=unique(core$pop.char))%do%{
    print(p)
    #calculate actual change and logit -transformed change in allele frequency. spring-fall so will be positive when spring freqeuncy is higher (same sign as gwas)
    delta=core[pop.char==p,.(diff.logit=logit.rf[season=="spring"]-logit.rf[season=="fall"], diff=rf[season=="spring"]-rf[season=="fall"],population=p), .( chr, pos )]
    return(delta)
}
names(deltas)=unique(core$pop.char)
save(deltas, file= "/mnt/pricey_2/priscilla/core20delta.rdat")

#move to rivanna and code there now