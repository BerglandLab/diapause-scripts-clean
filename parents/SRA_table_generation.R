library(data.table)
library(stringr)
sra<-fread("~/Box Sync/scripts-clean/parents/sraData.txt", header=F)
names(sra)=c("messy", "SRA")
sra[,id:=tstrsplit(messy, split="/")[[2]]]
sra[,loc:=tstrsplit(messy, split="/")[[1]]]

sra[loc=="12LN6", loc2:="12LN6_"]
sra[loc=="Ithica", loc2:=""]
sra[loc=="DGRP", loc2:="DGRP_"]
sra[loc=="12LN10", loc2:="12LN10_"]
sra[loc=="12BME10", loc2:="12BME10_"]

sra[,strain:=ifelse(is.na(loc2), id, paste0(loc2, id))]
geo<-fread("~/Box Sync/HSparents/strain_geography.csv")
a.lines<-fread("~/Box Sync/hybridSwarm/swarmA_lines.csv", header=F)[,V1]
b.lines<-fread("~/Box Sync/hybridSwarm/swarmB_lines.csv", header=F)[,V1]

sra$strain%in%geo$strain
sra$strain[45:50]

a<-merge(sra, geo, by="strain")
a[,messy:=NULL]
a[,id:=NULL]
a[,loc:=NULL]
a[,loc2:=NULL]
setnames(a, "actual location", "location")

a.sum<-a[, .(SRAs=paste(as.character(list(SRA)), sep=", ")), .(strain, geography, latitude, longitude, location)]

a.sum[, SRAs:=str_remove(SRAs, "c")]
a.sum[, SRAs:=str_remove_all(SRAs, "[:punct:]")]
a.sum[, SRAs:=str_remove(SRAs, "\)")]



a.sum[strain%in%a.lines, population:="A"]
a.sum[strain%in%b.lines, population:="B"]


write.table(a.sum[order(geography)][order(population)], "~/Box Sync/scripts-clean/parents/all_line_data_collapsed.txt", sep="\t", quote=F, row.names = F)
