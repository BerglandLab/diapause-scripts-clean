#get wolbachia read counts for all hybrid swarm individuals

ls *.bam > files.txt
while read file; do
  echo $file
  samtools idxstats $file > $file.stats
done<files.txt
