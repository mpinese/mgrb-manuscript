#!/bin/bash
set -euo pipefail

INFILE="${1}"
OUTFILE="${2}"

../../software/ensembl-vep/vep -v \
  -i "${INFILE}" --format vcf \
  --tab -o "${OUTFILE}" --force --no_stats \
  --everything --hgvsg --pick --no_escape \
  --offline --dir_cache ../../software/.vep/ --fasta ../../../resources/hs37d5x/reference_genome/hs37d5x.fa
