
library(data.table)
library(foreach)
library(rehh)
library(doMC)
registerDoMC(20)

#read ihs. data
a.raw<-fread("/scratch/pae3g/revisions/evolution/ihs.txt")

#read file info to loop through
files<-fread("/scratch/pae3g/genome-reconstruction/universal_input.txt")


    
    j<-foreach(pop=files$V1[1:3000], phenotype=files$V2[1:3000],draw=files$V3[1:3000], perm=files$V4[1:3000], model=files$V6[1:3000], .errorhandling="remove")%dopar%{
      print(paste(draw, perm, sep=","))
      #read gwas adn do some fixes
      load(paste("/scratch/pae3g/revisions/genesis_", phenotype, "_draw", draw, "_perm", perm, "_pop", pop, "_" , model, "_allsnpgrm_wolbachia_dropmissing.Rdat", sep=""))
      gwas<-assoc.results   
      gwas[,maf:=pmin(freq, 1-freq)]
      gwas<-gwas[maf>=0.05]
      gwas[,q:=log10(frank(Score.pval)/(length(Score.pval)+1))]
      setnames(gwas, c("chr", "pos"), c("CHR", "POSITION"))
      
      #merge with ihs data
      a<-merge(a.raw, gwas, by=c("CHR", "POSITION"), all.x=T)
      #currently coded so that ancestral allele is reference (0) and derived allele is alternate (1)
      #we want "derived" to be the pro-diapause allele
      #positive score means that ref is pro-diapause (aka ancestral is pro diapause)
      #so, when gwas score is POSTIVE, switch ihh_a and ihh_d before computing ihs
      
      a[Score>=0,ihh.a:=iHH_D]
      a[is.na(Score)|Score<0,ihh.a:=iHH_A]
      a[Score>=0,ihh.d:=iHH_A]
      a[is.na(Score)|Score<0,ihh.d:=iHH_D]
      
      #flip allele frequency if ancestral and derived are flipped
      a[Score>=0, freq_A:=1-freq_A]
      
      #delete original calls and rename
      a[,iHH_A:=NULL]
      a[,iHH_D:=NULL]
      setnames(a, c("ihh.a", "ihh.d"), c("iHH_A", "iHH_D"))
      a<-a[,.(CHR, POSITION, freq_A,  iHH_A, iHH_D, iES_Tang_et_al_2007, iES_Sabeti_et_al_2007)]
      
      b<-ihh2ihs(a) #note the default is discard maf < 0.05
      c<-as.data.table(b$ihs)
      
      d<-merge(gwas, c, by=c("CHR", "POSITION"))
      setnames(d, c("CHR","POSITION", "LOGPVALUE"), c("chr", "pos",'ihs.p'))
      
      x<-foreach (top = seq(from=-5, to=0, by=1)) %do% { 
        g.sum<-d[q>top,.(min.ihs=min(IHS, na.rm=T),
                        
                          med.ihs=median(IHS, na.rm=T), 
                          
                          max.ihs=max(IHS, na.rm=T),
                          n=.N)]
        g.sum[,perm:=perm]
        g.sum[,draw:=draw]
        g.sum[,top:=top]
        g.sum[,pheno:=phenotype]
        g.sum[,pop:=pop]
        g.sum[,model:=model]
        return(g.sum)
      }
      return(rbindlist(x))
    }
    j<-rbindlist(j)
    
    j[,test:="DGRP"]
    #save output   
    output=paste0("/scratch/pae3g/revisions/evolution/ihs_dgrp_universal_threshold_dropmissing.Rdat")
    print(output)
    save(j, file=output)
    
a.raw<-fread("/scratch/pae3g/revisions/evolution/BME_LN_LNPA.ihs.txt")



k<-foreach(pop=files$V1[1:3000], phenotype=files$V2[1:3000],draw=files$V3[1:3000], perm=files$V4[1:3000], model=files$V6[1:3000], .errorhandling="remove")%dopar%{
  print(paste(draw, perm, sep=","))
  #read gwas adn do some fixes
  load(paste("/scratch/pae3g/revisions/genesis_", phenotype, "_draw", draw, "_perm", perm, "_pop", pop, "_" , model, "_allsnpgrm_wolbachia_dropmissing.Rdat", sep=""))
  gwas<-assoc.results   
  gwas[,maf:=pmin(freq, 1-freq)]
  gwas<-gwas[maf>=0.05]
  gwas[,q:=log10(frank(Score.pval)/(length(Score.pval)+1))]
  setnames(gwas, c("chr", "pos"), c("CHR", "POSITION"))
  
  #merge with ihs data
  a<-merge(a.raw, gwas, by=c("CHR", "POSITION"), all.x=T)
  #currently coded so that ancestral allele is reference (0) and derived allele is alternate (1)
  #we want "derived" to be the pro-diapause allele
  #positive score means that ref is pro-diapause (aka ancestral is pro diapause)
  #so, when gwas score is POSTIVE, switch ihh_a and ihh_d before computing ihs
  
  a[Score>=0,ihh.a:=iHH_D]
  a[is.na(Score)|Score<0,ihh.a:=iHH_A]
  a[Score>=0,ihh.d:=iHH_A]
  a[is.na(Score)|Score<0,ihh.d:=iHH_D]
  
  #flip allele frequency if ancestral and derived are flipped
  a[Score>=0, freq_A:=1-freq_A]
  
  #delete original calls and rename
  a[,iHH_A:=NULL]
  a[,iHH_D:=NULL]
  setnames(a, c("ihh.a", "ihh.d"), c("iHH_A", "iHH_D"))
  a<-a[,.(CHR, POSITION, freq_A,  iHH_A, iHH_D, iES_Tang_et_al_2007, iES_Sabeti_et_al_2007)]
  
  b<-ihh2ihs(a) #note the default is discard maf < 0.05
  c<-as.data.table(b$ihs)
  
  d<-merge(gwas, c, by=c("CHR", "POSITION"))
  setnames(d, c("CHR","POSITION", "LOGPVALUE"), c("chr", "pos",'ihs.p'))
  
  x<-foreach (top = seq(from=-5, to=0, by=1)) %do% { 
    g.sum<-d[q>top,.(min.ihs=min(IHS, na.rm=T),
                      
                      med.ihs=median(IHS, na.rm=T), 
                      
                      max.ihs=max(IHS, na.rm=T),
                      n=.N)]
    g.sum[,perm:=perm]
    g.sum[,draw:=draw]
    g.sum[,top:=top]
    g.sum[,pheno:=phenotype]
    g.sum[,pop:=pop]
    g.sum[,model:=model]
    return(g.sum)
  }
  return(rbindlist(x))
}
k<-rbindlist(k)
k[,test:="North"]


#save output   
output=paste0("/scratch/pae3g/revisions/evolution/ihs_north_universal_threshold_dropmissing_bottomquantiles.Rdat")
print(output)
save(k, file=output)



#make plots
# load("/scratch/pae3g/revisions/evolution/ihs_dgrp_universal_threshold.Rdat")
# load("/scratch/pae3g/revisions/evolution/ihs_north_universal_threshold.Rdat")
#     
# j[,test:="DGRP"]
# k[,test:="Northern"]
# 
# i<-rbind(j,k)
# 
# i[perm!=0&pop=="both", group:="both-permuted"]
# i[perm==0&pop=="both", group:="both-observed"]
# i[perm==0&pop=="A", group:="A-observed"]
# i[perm==0&pop=="B", group:="B-observed"]
# 
# i[perm!=0&pop=="A", group:="A-permuted"]
# i[perm!=0&pop=="B", group:="B-permuted"]
# 
# i[,pheno2:=ifelse(pheno=="diapause.bin", "st. 8", "st. 10")]
# 
# i[top==-5, th:="Top 0.001%"]
# i[top==-4, th:="Top 0.01%"]
# i[top==-3, th:="Top 0.1%"]
# i[top==-2, th:="Top 1%"]
# i[top==-1, th:="Top 10%"]
# i[top==0, th:="all SNPs"]
# 
# i.melt<-melt(i, id.vars=c("group", "pheno", "pheno2", "th", "test", "pop", "perm", "draw"), measure.vars=c("min.ihs", "med.ihs", "max.ihs"))
# 
# y.sum<-i.melt[,.(med=median(value), q.025=quantile(value, .025), q.975=quantile(value, .975)), .(group, pheno, pheno2, th, test, pop, variable)]
# 
# 
# y.sum[,facet:=paste(pheno2, test, sep=": ")]
# y.sum[,facet:=factor(y.sum$facet, levels=c("st. 8: DGRP", "st. 10: DGRP", "st. 8: Northern", "st. 10: Northern"))]
# y.sum[,th:=factor(y.sum$th, levels=c("Top 0.001%", "Top 0.01%", "Top 0.1%", 'Top 1%', "Top 10%", "all SNPs"))]
# 
# a.plot<-ggplot(y.sum)+
#   geom_point(data=y.sum[variable=="med.ihs"], aes(x=pop, y=med, color=group), position=position_dodge(width=0.3))+
#   geom_errorbar(data=y.sum[variable=="med.ihs"], aes(x=pop, ymax=q.975, ymin=q.025, color=group), width=0.1, position=position_dodge(width=0.3))+
#   labs(x="", y="median IHS", color="")+
#   scale_color_manual(values=c("#39568CFF", "grey20", "#440154FF", "grey50", "lightseagreen", "grey80"))+
#   theme(legend.position = "none")+
#   facet_grid(th~facet, scales ="free")+
#   theme(axis.text.x=element_text(angle=45,hjust=1))
# 
# pdf("/scratch/pae3g/revisions/figures/ihs_threshold.pdf", height=8, width=8)
# a.plot
# dev.off()
