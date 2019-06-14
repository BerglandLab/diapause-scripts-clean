while read perm fdr top pthresh; do
plink1.9 \
--noweb \
--bfile /mnt/pricey_2/priscilla/final2  \
--clump /mnt/sammas_storage/bergland-lab/Priscilla/perm_summaries/gwas_p_socre_inv_id_perm"$perm."txt \
--clump-p1 $pthresh \
--clump-p2 0.1 \
--clump-kb 50 \
--clump-r2 0.2 \
--clump-field gwas.p \
--out /mnt/sammas_storage/bergland-lab/Priscilla/perm_summaries/perm$perm_top$top_ld0.4_50kb \
--allow-extra-chr
done < /mnt/pricey_2/priscilla/clump_p_thresholds_perms.txt
