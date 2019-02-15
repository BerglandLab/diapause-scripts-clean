
stems=paste("/scratch/pae3g/genome-reconstruction/final2_draw", c(1:100), "_replaced", sep="")
seeds=sample(1:100000, 100)
draw=c(1:100)
perm=c(101:200)

input<-data.table(fn=rep(stems, times=100, each=1), pop="both", phenotype="diapause.bin9", draws=rep(draw, times=100, each=1), perms=rep(perm, times=1, each=100), seed=rep(seeds, times=1, each=100))


write.table(input, "/scratch/pae3g/genome-reconstruction/seeded_permutation_input.txt", row.names=F, col.names=F, quote=F, sep="\t")

#adaptive permutation input

seeds=sample(1:1000000, 10000)
perm=c(1001:11000)

input<-data.table(pop="both", pheno="diapause.bin9", perm=perm, seed=seeds)

write.table(input, "/scratch/pae3g/genome-reconstruction/adaptive_permutation_input.txt", row.names=F, col.names=F, quote=F, sep="\t")

#stage 7 diapause input


stems=paste("/scratch/pae3g/genome-reconstruction/final2_draw", c(1:100), "_replaced", sep="")
seeds=sample(1:100000, 100)
draw=c(1:100)
perm=c(101:200)

input<-data.table(fn=rep(stems, times=100, each=1), pop="both", phenotype="diapause.bin", draws=rep(draw, times=100, each=1), perms=rep(perm, times=1, each=100), seed=rep(seeds, times=1, each=100))


write.table(input, "/scratch/pae3g/genome-reconstruction/seeded_permutation_input2.txt", row.names=F, col.names=F, quote=F, sep="\t")

#make input for actual and for last line that was #10,000 and couldn't run

stems=paste("/scratch/pae3g/genome-reconstruction/final2_draw", c(1:100), "_replaced", sep="")
seeds=sample(1:100000, 100)
draw=c(1:100)
perm=c(101:200)

input<-data.table(fn=rep(stems, times=100, each=1), pop="both", phenotype="diapause.bin", draws=rep(draw, times=100, each=1), perms=rep(perm, times=1, each=100), seed=rep(seeds, times=1, each=100))

input<-input[10000,]
input2<-data.table(fn=stems, pop="both", phenotype="diapause.bin", draws=c(1:100), perms=rep(0, times=100), seed=rep(1, times=100))

input2=rbind(input, input2)

write.table(input2, "/scratch/pae3g/genome-reconstruction/seeded_permutation_input3.txt", row.names=F, col.names=F, quote=F, sep="\t")

