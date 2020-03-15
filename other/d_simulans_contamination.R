# d simulans in 2019
library(data.table)
library(lubridate)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
w<-fread("~/Box Sync/Job applications/2019/alyssa_zap_counts.csv")
w[,coll.date:=mdy(`Date Collected`)]

w.sum<-w[,.(prop.dsim=sum(COUNT[species=="D. simulans" | species=="Drosophila simulans"], na.rm=T)/sum(COUNT[species=="D. simulans" | species=="Drosophila simulans"|species=="D. melanogaster"], na.rm=T), total=sum(COUNT[species=="D. simulans" | species=="Drosophila simulans"|species=="D. melanogaster"], na.rm=T)), .(coll.date)]

w.sum[,year:=year(coll.date)]
w.sum[,j:=yday(coll.date)]

a<-ggplot(data=w.sum, aes(x=coll.date, y=prop.dsim))+geom_point()+geom_line()+scale_x_date(limits=as.Date(c('2017-06-15', '2017-12-15')))+labs(x=NULL, y="Proportion D. sim", title="2017")
b<-ggplot(data=w.sum, aes(x=coll.date, y=prop.dsim))+geom_point()+geom_line()+scale_x_date(limits=as.Date(c('2018-06-15', '2018-12-15')))+labs(x=NULL, y="Proportion D. sim", title="2018")

plot_grid(a,b, nrow=2)
