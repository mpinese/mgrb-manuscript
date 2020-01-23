#!/bin/bash
set -euo pipefail

REFERENCE="../../../resources/hs37d5x/reference_genome/hs37d5x.fa"        # Local
#REFERENCE="../../../../../resources/hs37d5x/reference_genome/hs37d5x.fa"  # Raijin

mkdir -p tmp

gzip -dc CosmicCodingMuts.vcf.gz | \
  ./bin/vt normalize -m -r "${REFERENCE}" - | \
  awk 'BEGIN {FS="\t"; OFS="\t"; print "chrom", "pos", "ref", "alt", "cosmic.id", "cosmic.count"}  /^[^#]/ {cnt=$8; gsub(/.*CNT=/, "", cnt); gsub(/;.*/, "", cnt); print $1, $2, $4, $5, $3, cnt}' \
  > ../04_cosmic_83.tsv

