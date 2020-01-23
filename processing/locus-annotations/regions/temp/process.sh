gzip -dc ../bad.*.bed.gz | awk 'BEGIN{OFS="\t"} {if ($2<5) {$2=0} else {$2-=5}; print $1, $2, $3+5}' | cat /dev/stdin par.bed noncanonical.bed | sort -k1,1 -k2,2n | bedtools merge > tier3.bed
bedtools subtract -a canonical.bed -b tier3.bed > tier12.bed
bedtools intersect -a tier12.bed -b ../good.HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed.gz > tier1.bed
bedtools subtract -a tier12.bed -b tier1.bed > tier2.bed

awk '{total += $3-$2} END {print total/1e6}' < tier1.bed
awk '{total += $3-$2} END {print total/1e6}' < tier2.bed
awk '{total += $3-$2} END {print total/1e6}' < tier3.bed
awk '{total += $3-$2} END {print total/1e6}' < genome.bed

bedtools intersect -a tier1.bed -b canonical.bed | awk '{total += $3-$2} END {print total/1e6}'
bedtools intersect -a tier2.bed -b canonical.bed | awk '{total += $3-$2} END {print total/1e6}'
bedtools intersect -a tier3.bed -b canonical.bed | awk '{total += $3-$2} END {print total/1e6}'
awk '{total += $3-$2} END {print total/1e6}' < canonical.bed

bedtools intersect -a tier1.bed -b ../../code/ccds-20171121.bed | awk '{total += $3-$2} END {print total/1e6}'
bedtools intersect -a tier2.bed -b ../../code/ccds-20171121.bed | awk '{total += $3-$2} END {print total/1e6}'
bedtools intersect -a tier3.bed -b ../../code/ccds-20171121.bed | awk '{total += $3-$2} END {print total/1e6}'
awk '{total += $3-$2} END {print total/1e6}' < ../../code/ccds-20171121.bed
