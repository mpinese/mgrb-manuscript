#!/bin/bash
set -euo pipefail

mkdir -p ../converted_data

# Create a BED of coding regions.  CADD scores will be subset to these,
# otherwise I've had a lot of difficulty getting the whole genome into
# Hail.  Noncoding CADD scores are of little interest anyway.
# Use a window of +/- 10 kb around the UCSC Ensembl coding regions for
# safety.
awk 'BEGIN{FS="\t";OFS="\t"} {print $1, 0, $2}' < ../../../../resources/hs37d5x/reference_genome/hs37d5x.fa.fai > 02_genome.bed
mysql -NB -h genome-mysql.cse.ucsc.edu hg19 -u genome -A -e 'select chrom, txStart, txEnd from ensGene' | \
  awk 'BEGIN {FS="\t";OFS="\t"} {chrom=$1; sub(/^chr/, "", chrom); if (chrom == "M") { chrom = "MT" }; $2 -= 10000; $3 += 10000; if ($2 < 0) { $2 = 0 }; print chrom, $2, $3}' | \
  sort -k1,1 -k2,2n | \
  bedtools intersect -a - -b 02_genome.bed | \
  bedtools merge | bgzip > 02_target_regions.bed.gz
tabix 02_target_regions.bed.gz

# Subset the downloaded scores to these coding regions
#parallel ./tab2vcf2.sh {} ../converted_data/{/.}.cds.vcf.bgz 02_target_regions.bed.gz ::: ../source_data/1000G_phase3.tsv.gz ../source_data/ESP6500SI.tsv.gz ../source_data/ExAC_r0.3.tsv.gz
parallel ./tab2vcf.sh {1} ../converted_data/{1/.}.cds.vcf.bgz 02_target_regions.bed.gz {2} ::: ../source_data/whole_genome_SNVs.tsv.gz ::: X $(seq 1 22) Y MT
