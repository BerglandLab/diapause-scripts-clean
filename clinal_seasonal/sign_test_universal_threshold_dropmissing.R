library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(foreach)
library(doMC)
registerDoMC(16)
library(viridis)
library(ggbeeswarm)
library(Rmisc)

print("bergland 2014")

load("/scratch/pae3g/oldscratch_recovered/evolution/6d_data.Rdata")
b<-as.data.table(p)
b[,maf:=pmin(f.hat, 1-f.hat)]
b<-b[maf>=0.05]
b<-b[clinal.beta<3&clinal.beta>(-3)]
b<-b[sfsfsfX.beta<3&sfsfsfX.beta>(-3)]
b[,clinal.q:=log10(frank(clinal.p)/(length(clinal.p)+1))]
b[,sfsfsfX.q:=log10(frank(sfsfsfX.p)/(length(sfsfsfX.p)+1))]

files<-fread("/scratch/pae3g/genome-reconstruction/universal_input.txt")



#read in files
y<-foreach(pop=files$V1[1:3000], phenotype=files$V2[1:3000],draw=files$V3[1:3000], perm=files$V4[1:3000], model=files$V6[1:3000])%dopar%{
    print(paste(draw, perm, sep=","))
    #read gwas adn do some fixes
    load(paste("/scratch/pae3g/revisions/genesis_", phenotype, "_draw", draw, "_perm", perm, "_pop", pop, "_" , model, "_allsnpgrm_wolbachia_dropmissing.Rdat", sep=""))
    gwas<-assoc.results
    gwas[,maf:=pmin(freq, 1-freq)]
    gwas<-gwas[maf>=0.05]
    gwas[,q:=log10(frank(Score.pval)/(length(Score.pval)+1))]
    gwas<-merge(gwas, b[,.(chr, pos, clinal.beta, clinal.p,clinal.q,  maf, sfsfsfX.p, sfsfsfX.beta, sfsfsfX.q)], by=c("chr", "pos"))
    gwas[,TT.cline:=ifelse(sign(Score)==1 & sign(clinal.beta)==(-1), T, F)]
    gwas[,TF.cline:=ifelse(sign(Score)==1 & sign(clinal.beta)==1, T, F)]
    gwas[,FT.cline:=ifelse(sign(Score)==(-1) & sign(clinal.beta)==(-1), T, F)]
    gwas[,FF.cline:=ifelse(sign(Score)==(-1) & sign(clinal.beta)==1, T, F)]

    gwas[,TT.seas:=ifelse(sign(Score)==1 & sign(sfsfsfX.beta)==1, T, F)]
    gwas[,TF.seas:=ifelse(sign(Score)==1 & sign(sfsfsfX.beta)==(-1), T, F)]
    gwas[,FT.seas:=ifelse(sign(Score)==(-1) & sign(sfsfsfX.beta)==1, T, F)]
    gwas[,FF.seas:=ifelse(sign(Score)==(-1) & sign(sfsfsfX.beta)==(-1), T, F)]

    #bergland cline is based on alternate alleles and needs to be flipped. 
    gwas[,ps.cline:=-1*Score.Stat*clinal.beta]
    #bergland positive clinal beta means alternate allele is higher in the fall (same direction as gwas)
    gwas[,ps.seas:=1*Score.Stat*sfsfsfX.beta]


    cline.top<-foreach (top = seq(from=-5, to=0, by=1)) %do% {
        a<-gwas[q<top,.(
            n=.N,
            TT=sum(TT.cline, na.rm=T),
            TF=sum(TF.cline, na.rm=T),
            FT=sum(FT.cline, na.rm=T),
            FF=sum(FF.cline, na.rm=T),
            poly=sum(ps.cline)
        )]
        a[,top:=top]
    }
    cline.top<-rbindlist(cline.top)
    cline.top[,or:=as.numeric(TT)*as.numeric(FF)/(as.numeric(FT)*as.numeric(TF))]
    cline.top[, perm:=perm]
    cline.top[,draw:=draw]
    cline.top[,test:="clinal"]
    cline.top[,pheno:=phenotype]
    cline.top[,model:=model]
    cline.top[,pop:=pop]

    

    seas.top<-foreach (top = seq(from=-5, to=0, by=1)) %do% {
        a<-gwas[q<top,.(
            n=.N,
            TT=sum(TT.seas, na.rm=T),
            TF=sum(TF.seas, na.rm=T),
            FT=sum(FT.seas, na.rm=T),
            FF=sum(FF.seas, na.rm=T),
            poly=sum(ps.seas)
        )]
        a[,top:=top]
    }

    seas.top<-rbindlist(seas.top)
    seas.top[,or:=as.numeric(TT)*as.numeric(FF)/(as.numeric(FT)*as.numeric(TF))]
    seas.top[, perm:=perm]
    seas.top[,draw:=draw]
    seas.top[,test:="seasonal"]
    seas.top[,pheno:=phenotype]
    seas.top[,model:=model]
    seas.top[,pop:=pop]
    

    return(rbind(cline.top, seas.top))
}

y<-rbindlist(y)

save(y, file="/scratch/pae3g/revisions/evolution/bergland2014_sign_universal_threshold_dropmissing.Rdata")

print("machado 2019")

cline<-fread("/scratch/pae3g/revisions/evolution/east_coast_cline_V2_clean.txt")
seas<-fread("/scratch/pae3g/revisions/evolution/seas_glm_NOTswitch_clean.txt")   
freqs<-fread("/scratch/pae3g/oldscratch_recovered/evolution/east_coast_cline_V2_allele_freqs.txt")

x<-merge(cline, seas, by=c("chr", "pos"))
x<-merge(x, freqs, by=c("chr", "pos"))

x<-x[f.hat>=0.05&f.hat<=0.95]
x<-x[clinal.beta<3&clinal.beta>(-3)]
x<-x[seas.beta<3&seas.beta>(-3)]





y<-foreach(pop=files$V1[1:3000], phenotype=files$V2[1:3000],draw=files$V3[1:3000], perm=files$V4[1:3000], model=files$V6[1:3000], .errorhandling="remove")%dopar%{
    print(paste(draw, perm, sep=","))
    #read gwas adn do some fixes
    load(paste("/scratch/pae3g/revisions/genesis_", phenotype, "_draw", draw, "_perm", perm, "_pop", pop, "_" , model, "_allsnpgrm_wolbachia_dropmissing.Rdat", sep=""))
    gwas<-assoc.results
    gwas[,maf:=pmin(freq, 1-freq)]
    gwas<-gwas[maf>=0.05]
    gwas[,q:=log10(frank(Score.pval)/(length(Score.pval)+1))]
    gwas<-merge(gwas, x[,.(chr, pos, clinal.beta, clinal.p, seas.beta)], by=c("chr", "pos"))
    gwas[,TT.cline:=ifelse(sign(Score)==1 & sign(clinal.beta)==(-1), T, F)]
    gwas[,TF.cline:=ifelse(sign(Score)==1 & sign(clinal.beta)==1, T, F)]
    gwas[,FT.cline:=ifelse(sign(Score)==(-1) & sign(clinal.beta)==(-1), T, F)]
    gwas[,FF.cline:=ifelse(sign(Score)==(-1) & sign(clinal.beta)==1, T, F)]
    
    gwas[,TT.seas:=ifelse(sign(Score)==1 & sign(seas.beta)==(-1), T, F)]
    gwas[,TF.seas:=ifelse(sign(Score)==1 & sign(seas.beta)==(1), T, F)]
    gwas[,FT.seas:=ifelse(sign(Score)==(-1) & sign(seas.beta)==(-1), T, F)]
    gwas[,FF.seas:=ifelse(sign(Score)==(-1) & sign(seas.beta)==(1), T, F)]
    
    gwas[,ps.cline:=-1*Score.Stat*clinal.beta]
    #note that the seasonal polygenic score should have had a *-1 in it and needs to be flipped for plotting
    gwas[,ps.seas:=1*Score.Stat*seas.beta]
    
    
    cline.top<-foreach (top = seq(from=-5, to=0, by=1)) %do% {
        a<-gwas[q<top,.(
            n=.N,
            TT=sum(TT.cline, na.rm=T),
            TF=sum(TF.cline, na.rm=T),
            FT=sum(FT.cline, na.rm=T),
            FF=sum(FF.cline, na.rm=T),
            poly=sum(ps.cline)
        )]
        a[,top:=top]
    }
    
    cline.top<-rbindlist(cline.top)
    cline.top[,or:=as.numeric(TT)*as.numeric(FF)/(as.numeric(FT)*as.numeric(TF))]
    cline.top[, perm:=perm]
    cline.top[,draw:=draw]
    cline.top[,test:="clinal"]
    cline.top[,pheno:=phenotype]
    cline.top[,model:=model]
    cline.top[,pop:=pop]
    
    seas.top<-foreach (top = seq(from=-5, to=0, by=1)) %do% {
        a<-gwas[q<top,.(
            n=.N,
            TT=sum(TT.seas, na.rm=T),
            TF=sum(TF.seas, na.rm=T),
            FT=sum(FT.seas, na.rm=T),
            FF=sum(FF.seas, na.rm=T),
            poly=sum(ps.seas)
        )]
        a[,top:=top]
    }
    
    seas.top<-rbindlist(seas.top)
    seas.top[,or:=as.numeric(TT)*as.numeric(FF)/(as.numeric(FT)*as.numeric(TF))]
    seas.top[, perm:=perm]
    seas.top[,draw:=draw]
    seas.top[,test:="seasonal"]
    seas.top[,pheno:=phenotype]
    seas.top[,model:=model]
    seas.top[,pop:=pop]
    
    
    return(rbind(cline.top, seas.top))
}

y<-rbindlist(y)

save(y, file="/scratch/pae3g/revisions/evolution/bergland2019_sign_universal_threshold_dropmissing.Rdata")


print("individual populations")
load("/scratch/pae3g/oldscratch_recovered/evolution/core20delta.rdat")
pops<-names(deltas)




y<-foreach(pop=files$V1[1:3000], phenotype=files$V2[1:3000],draw=files$V3[1:3000], perm=files$V4[1:3000], model=files$V6[1:3000], .errorhandling="remove")%dopar%{
    print(paste(draw, perm, sep=","))
    #read gwas adn do some fixes
    load(paste("/scratch/pae3g/revisions/genesis_", phenotype, "_draw", draw, "_perm", perm, "_pop", pop, "_" , model, "_allsnpgrm_wolbachia_dropmissing.Rdat", sep=""))
    gwas<-assoc.results
    gwas[,maf:=pmin(freq, 1-freq)]
    gwas<-gwas[maf>=0.05]
    gwas[,q:=log10(frank(Score.pval)/(length(Score.pval)+1))]

    pop.test<-foreach(p=pops)%do%{
        print(p)
        seas.top<-foreach (top = c(-5:0)) %do% {
            
            l<-merge(gwas[q<=top], deltas[[p]][is.finite(diff.logit)], by=c("chr", "pos"))
            l[,TT:=ifelse(sign(Score.Stat)==1 & sign(diff.logit)==(1), T, F)]
            l[,TF:=ifelse(sign(Score.Stat)==(-1) & sign(diff.logit)==1, T, F)]
            l[,FT:=ifelse(sign(Score.Stat)==(1) & sign(diff.logit)==(-1), T, F)]
            l[,FF:=ifelse(sign(Score.Stat)==(-1) & sign(diff.logit)==(-1), T, F)]
            l[,ps:=Score.Stat*diff.logit]
            l[,top:=top]
            return(l[,.(or=sum(as.numeric(TT), na.rm=T)*sum(as.numeric(FF), na.rm=T)/(sum(as.numeric(TF), na.rm=T)*sum(as.numeric(FT), na.rm=T)), 
                        poly=sum(ps)), .(population, top)])
            
        }
        return(rbindlist(seas.top))
    }
    
    
    
    pop.test<-rbindlist(pop.test)
    pop.test[,pheno:=phenotype]
    pop.test[,pop:=pop]
    pop.test[, perm:=perm]
    pop.test[,draw:=draw]
    pop.test[,model:=model]
    return(pop.test)
}

y<-rbindlist(y)

save(y, file="/scratch/pae3g/revisions/evolution/single_population_sign_universal_threshold_dropmissing.Rdata")




