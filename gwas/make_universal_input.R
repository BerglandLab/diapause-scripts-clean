#make universal input file


i<-fread("/scratch/pae3g/genome-reconstruction/1000perm_100imp_stage7.txt")
i[,V1:=NULL]


i2=copy(i)
i2[,V3:="diapause.bin9"]

i3<-rbind(i[1:200], i2[1:200])
i3[,V2:="A"]
i4<-copy(i3)
i4[,V2:="B"]
i5<-rbind(i, i2, i3, i4)

i5[,V7:="nonloco"]
i6<-copy(i5)
i6[,V7:="loco"]

i7<-rbind(i5, i6)

write.table(i7, "/scratch/pae3g/genome-reconstruction/universal_input.txt", quote=F, row.names=F, col.names=F, sep="\t")
