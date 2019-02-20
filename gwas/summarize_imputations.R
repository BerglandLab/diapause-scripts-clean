library(data.table)
library(foreach)
library(cowplot)
library(doMC)
registerDoMC(16)

foreach(perm=c(101:200))%do%{
    print(perm)
    x<-foreach(draw=c(1:100))%dopar%{
        print(draw)
        return(fread(paste("/scratch/pae3g/final_reconstruction2/genesis_diapause.bin9~generation+temp.rack.cal+photoperiod_both_draw", draw, "_perm", perm, "_replaced.txt", sep="")))
    }
    x<-rbindlist(x)
    y<-x[,.(mean.p=mean(Score.pval, na.rm=T), med.p=median(Score.pval, na.rm=T), min.p=min(Score.pval, na.rm=T), sd.pval=sd(Score.pval, na.rm=T), mean.Score=mean(Score, na.rm=T), med.Score=median(Score, na.rm=T), sd.Score=sd(Score, na.rm=T), mean.Score.Stat=mean(Score.Stat, na.rm=T), med.Score.Stat=median(Score.Stat, na.rm=T), sd.Score.Stat=sd(Score.Stat, na.rm=T)), .(snpID, chr, pos)]
    write.table(y, paste("/scratch/pae3g/final_reconstruction2/perm", perm,"_100_draw_summary.txt", sep=""), sep="\t", row.names=F, col.names=T, quote=F)
}

