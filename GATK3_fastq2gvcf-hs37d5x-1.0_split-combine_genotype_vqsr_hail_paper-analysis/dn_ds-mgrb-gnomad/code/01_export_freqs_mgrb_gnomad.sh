#!/bin/bash
set -euo pipefail

mkdir -p tmp
cp goi_list.txt tmp/gois.txt

awk 'BEGIN{FS="\t";OFS="\t"} {print $1, 0, $2}' < ../../../../resources/hs37d5x/reference_genome/hs37d5x.fa.fai > tmp/genome.bed

query="SELECT chrom, txStart, txEnd FROM refGene WHERE name2 IN ("
while read gene_symbol; do
  query="${query} '${gene_symbol}', "
done < <(tail -n+2 tmp/gois.txt | cut -f3 | sort | uniq)
query="${query%, } );"

mysql -NB -h genome-mysql.cse.ucsc.edu hg19 -u genome -A -e "${query}" | \
  awk 'BEGIN {FS="\t";OFS="\t"} {chrom=$1; sub(/^chr/, "", chrom); if (chrom == "M") { chrom = "MT" }; $2 -= 10000; $3 += 10000; if ($2 < 0) { $2 = 0 }; print chrom, $2, $3}' | \
  sort -k1,1 -k2,2n | \
  bedtools intersect -a - -b tmp/genome.bed | \
  bedtools merge | bgzip > tmp/gois.bed.gz
tabix tmp/gois.bed.gz


. 01h_download_gnomad_wgs.sh
./bin/pyhail.sh 01h_export_gois_freqs_gnomadwgs_all.py

./bin/pyhail.sh 01h_export_gois_freqs_mgrb_all.py
