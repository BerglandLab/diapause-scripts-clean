### libraries
	library("data.table")
	library("foreach")
		
### pull in SLURM ARRAY BATCH ID
	args.vec <- as.numeric(commandArgs(trailing=T)) + 1
#    args.vec <- 1
    
### load gwas file list
	gwas.fl <- data.table(fl=list.files("/nv/vol186/bergland-lab/Priscilla/perm_summaries", "adaptive_perms_perm", full.names=T))
	gwas.fl[grepl("perm0", fl), type:="orig"]
	gwas.fl[!grepl("perm0", fl), type:="perm"]
	
### load pvalue vector files
	### this is for the seasonal set + population level Fisher tests.
		#pval.fl <- data.table(fl=c(list.files("/scratch/aob2x/diapause_gwas", "fisher_exactJ", full.names=T),
		#							"/scratch/aob2x/diapause_gwas/mel_all_paired20_2sample_caF_popyear_4switch.f_s.glm"))
	
	### seasonal & clinal set
		pval.fl <- data.table(fl=c("/nv/vol186/bergland-lab/Priscilla/east_coast_cline_V2_clean.txt",
									"/scratch/aob2x/diapause_gwas/mel_all_paired20_2sample_caF_popyear_4switch.f_s.glm"))
									
### make full combo matrix
	combos <- as.data.table(expand.grid(pval.fl$fl, gwas.fl$fl)) 
	combos[,Var1:=as.character(Var1)]
	combos[,Var2:=as.character(Var2)]
	
### simplify combos
	#combos <- combos[grepl("paired20", Var1)]
	
### load in selected GWAS
	gwas.obs <- fread(combos[args.vec]$Var2, sep="\t")
	
	if(grepl("perm0", combos[args.vec]$Var2)) {
		gwas.obs <- gwas.obs[,c("V2","V3","V10"),with=F]
		setnames(gwas.obs, c("V2","V3","V10"), c("chr", "pos", "gwas.p"))
	} else {
		setnames(gwas.obs, c("pval.percentile"), c("gwas.p"))
	}

### load in selected seasonal/clinal set (note - the object is called seasonal set, but this is just a carry over; columns are being renamed "seas.p" jsut so that I don't have to rename below.
	seasonal.pval <- fread(combos[args.vec]$Var1)
	
	if(any(grepl("minp2", names(seasonal.pval)))) {
		setnames(seasonal.pval, c("chrom", "minp2"), c("chr", "seas.p"))
	} else if(any(grepl("minp", names(seasonal.pval)))) {
		setnames(seasonal.pval, c("chrom", "minp"), c("chr", "seas.p"))
	} else if(any(grepl("clinal.p", names(seasonal.pval)))) {
		setnames(seasonal.pval, "clinal.p", "seas.p") 
	} else if(any(names(seasonal.pval)=="chrom")) {
		setnames(seasonal.pval, c("chrom"), c("chr"))
	} 

### merge orignal GWAS & selected seasonal set
	setkey(gwas.obs, chr, pos)
	setkey(seasonal.pval, chr, pos)
	
	m.orig <- merge(gwas.obs, seasonal.pval)
	m.orig <- m.orig[gwas.p!=0][seas.p!=0]

### generate rank-normalized scores
	m.orig[,gwas.q := (rank(gwas.p)/(length(gwas.p)+1))]
	m.orig[,seas.q := (rank(seas.p)/(length(seas.p)+1))]
	
### convert to Z
	m.orig[,gwas.Zq := qnorm(gwas.q, 0, 1)]
	m.orig[,seas.Zq := qnorm(seas.q, 0, 1)]
	
	m.orig[,gwas.Zp := qnorm(gwas.p, 0, 1)]
	m.orig[,seas.Zp := qnorm(seas.p, 0, 1)]
	
### merge p-values
	m.orig[,st.Zp := (seas.Zp+gwas.Zp)/sqrt(2)]
	m.orig[,st.Zp.p:=pnorm(st.Zp, 0, 1)]
	m.orig[,st.Zp.q:=p.adjust(st.Zp.p, "fdr")]

	m.orig[,st.Zq := (seas.Zq+gwas.Zq)/sqrt(2)]
	m.orig[,st.Zq.p:=pnorm(st.Zq, 0, 1)]
	m.orig[,st.Zq.q:=p.adjust(st.Zq.p, "fdr")]

### save
	write.csv(m.orig, paste("/scratch/aob2x/diapause_gwas/mvn_output/stouffers.job_", args.vec, ".csv", sep=""),
				quote=F, row.names=F)
	
### summary
	m.orig.ag <- data.table(i=args.vec, 
				min.st.Zp.p=min(m.orig$st.Zp.p), min.st.Zp.q=min(m.orig$st.Zp.q),
				min.st.Zq.p=min(m.orig$st.Zq.p), min.st.Zq.q=min(m.orig$st.Zq.q))
				

	write.csv(m.orig.ag, 
				paste("/scratch/aob2x/diapause_gwas/mvn_output/stouffers.summ.job_", args.vec, ".csv", sep=""),
				quote=F, row.names=F)
	
	
	
	
	### library
	library(data.table)
	library(foreach)
	library(stringr)
	
	### file-list
	fl.dt <- data.table(fl=list.files("/nv/vol186/bergland-lab/diapause_gwas/mvn_output/", 
	                                  "stouffers.summ.job_", full.names=T))
	fl.dt[,i:= as.numeric(gsub(".csv", "", do.call("rbind", strsplit(unlist(lapply(strsplit(fl.dt$fl, "/"), function(x) rev(x)[1])), "_"))[,2]))]
	setkey(fl.dt, i)
	
	### load gwas file list
	gwas.fl <- data.table(fl=list.files("/nv/vol186/bergland-lab/Priscilla/perm_summaries", "adaptive_perms_perm", full.names=T))
	gwas.fl[grepl("perm0", fl), type:="orig"]
	gwas.fl[!grepl("perm0", fl), type:="perm"]
	
	### load pvalue vector files
	### this is for the seasonal set + population level Fisher tests.
	#pval.fl <- data.table(fl=c(list.files("/scratch/aob2x/diapause_gwas", "fisher_exactJ", full.names=T),
	#							"/scratch/aob2x/diapause_gwas/mel_all_paired20_2sample_caF_popyear_4switch.f_s.glm"))
	
	### seasonal & clinal set
	pval.fl <- data.table(fl=c("/nv/vol186/bergland-lab/Priscilla/east_coast_cline_V2_clean.txt",
	                           "/scratch/aob2x/diapause_gwas/mel_all_paired20_2sample_caF_popyear_4switch.f_s.glm"))
	
	### make full combo matrix
	combos <- as.data.table(expand.grid(pval.fl$fl, gwas.fl$fl))
	combos[,Var1:=as.character(Var1)]
	combos[,Var2:=as.character(Var2)]
	combos[,i:=1:dim(combos)[1]]
	setkey(combos ,i)
	
	### merge
	combos <- merge(combos, fl.dt, all.x=T)	
	combos[,seasonalSet := unlist(lapply(strsplit(combos$Var1, "_"), function(x) paste(gsub(".coef", "", x[c(4,5)]), collapse="_")))]
	combos[,gwasSet := as.numeric(str_match(str_match(combos$Var2, "perm[0-9]{1,}"), "[0-9]{1,}"))]
	
	### remove NAs
	combos <- combos[!is.na(fl)]
	
	### load
	o <- foreach(i=1:dim(combos)[1])%do%{
	    if(i%%25==0) print(paste(i, dim(combos)[1], sep=" / "))
	    temp <- fread(combos$fl[i])
	    cbind(temp, combos[i,c("seasonalSet", "gwasSet", "i"), with=F])
	}
	o <- rbindlist(o)
	
	### some basic summaries: basically, how many times did the minimum combined p-value in the observed data fall lower than the minimum across the permutations?
	o.ag <- o[,list(pr.min.st.Zp.p = mean(min.st.Zp.p[gwasSet==0] < min.st.Zp.p[gwasSet!=0], na.rm=T),
	                pr.min.st.Zp.q = mean(min.st.Zp.q[gwasSet==0] < min.st.Zp.q[gwasSet!=0], na.rm=T),
	                pr.min.st.Zq.p = mean(min.st.Zq.p[gwasSet==0] < min.st.Zq.p[gwasSet!=0], na.rm=T),
	                pr.min.st.Zq.q = mean(min.st.Zq.q[gwasSet==0] < min.st.Zq.q[gwasSet!=0], na.rm=T)), 
	          list(seasonalSet)]
	
	### critical values
	o.crit <- o[, list(quan.min.st.Zp.p = quantile(min.st.Zp.p[gwasSet!=0], .05, na.rm=T),
	                   min.st.Zp.p = min(min.st.Zp.p[gwasSet==0], na.rm=T),
	                   quan.min.st.Zq.p = quantile(min.st.Zq.p[gwasSet!=0], .05, na.rm=T),
	                   min.st.Zq.p = min(min.st.Zq.p[gwasSet==0], na.rm=T)), 
	            list(seasonalSet)]
	
	### load in sites passing critical values
	dat.sig <- foreach(ii=combos[gwasSet==0]$i, .errorhandling="remove")%do%{
	    print(ii)
	    dat <- fread(paste("/scratch/aob2x/diapause_gwas/mvn_output/stouffers.job_", ii, ".csv", sep=""))
	    dat[st.Zp.p <= o.crit[seasonalSet==combos[i==ii]$seasonalSet]$quan.min.st.Zp.p]
	}
	