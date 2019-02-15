
for swarm in A B; do
    #cd /scratch/pae3g/genome-reconstruction/output/cage$swarm
    #ls -d */ | cut -f1 -d'/' > /scratch/pae3g/genome-reconstruction/output/cage$swarm.txt
    while read i; do
	cd /scratch/pae3g/genome-reconstruction/output/cage$swarm/$i
	for x in *.tar.gz; do tar -xzf $x --wildcards --no-anchored '*.bam*'; done
	for num in {1..10}; do
	    cp $num.bam /scratch/pae3g/genome-reconstruction/output/bams$swarm/$i.$num.$swarm.bam
	    cp $num.bam.bai /scratch/pae3g/genome-reconstruction/output/bams$swarm/$i.$num.$swarm.bam.bai
	done
    done < /scratch/pae3g/genome-reconstruction/output/cage$swarm.txt
done

ls /scratch/pae3g/genome-reconstruction/output/bamsB/*.bam > /scratch/pae3g/genome-reconstruction/output/bamsB.list


ls /scratch/pae3g/genome-reconstruction/output/bamsA/*.bam > /scratch/pae3g/genome-reconstruction/output/bamsA.list

