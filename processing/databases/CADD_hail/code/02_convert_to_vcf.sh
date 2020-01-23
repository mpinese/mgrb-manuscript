#!/bin/bash
set -euo pipefail

mkdir -p ../converted_data

# Create a BED of coding regions.  CADD scores will be subset to these,
# otherwise I've had a lot of difficulty getting the whole genome into
# Hail.  Noncoding CADD scores are of little interest anyway.
# Use a window of +/- 10 kb around the UCSC Ensembl coding regions for
# safety.
awk 'BEGIN{FS="\t";OFS="\t"} {print $1, 0, $2}' < ../../../resources/hs37d5x/reference_genome/hs37d5x.fa.fai > 02_genome.bed
mysql -NB -h genome-mysql.cse.ucsc.edu hg19 -u genome -A -e 'select chrom, txStart, txEnd from ensGene' | \
  awk 'BEGIN {FS="\t";OFS="\t"} {chrom=$1; sub(/^chr/, "", chrom); if (chrom == "M") { chrom = "MT" }; $2 -= 10000; $3 += 10000; if ($2 < 0) { $2 = 0 }; print chrom, $2, $3}' | \
  sort -k1,1 -k2,2n | \
  bedtools intersect -a - -b 02_genome.bed | \
  bedtools merge | bgzip > 02_target_regions.bed.gz
tabix 02_target_regions.bed.gz

# Further split the coding regions into shards.  Grr.
#mkdir -p tmp
#gzip -dc 02_target_regions.bed.gz | \
#  awk 'BEGIN{FS="\t";OFS="\t";chrom="";outidx=0} {if (chrom!=$1 || size > 10000000) {outidx++; chrom=$1;size=0}; size+=$3-$2; print $0 > "tmp/02_target_regions_split"outidx".bed"}'
#
#parallel bgzip ::: tmp/*split*.bed
#parallel tabix ::: tmp/*split*.bed.gz

# Split by dataset and chromosome
parallel \
  --rpl '{/..} s:.*/::; s:\.[^/.]+$::; s:\.[^/.]+$::;' \
  -j14 \
  ./tab2vcf.sh {1} ../converted_data/{1/..}.{2}.vcf.bgz 02_target_regions.bed.gz {2} \
  ::: ../source_data/*.tsv.gz ::: X $(seq 1 22) Y MT

# Split by dataset and shard
#parallel \
#  --rpl '{/..} s:.*/::; s:\.[^/.]+$::; s:\.[^/.]+$::;' \
#  -j14 \
#  ./tab2vcf2.sh {1} ../converted_data/{1/..}.{2/.}.vcf.bgz {2} \
#  ::: ../source_data/*.tsv.gz ::: tmp/*split*.bed.gz
