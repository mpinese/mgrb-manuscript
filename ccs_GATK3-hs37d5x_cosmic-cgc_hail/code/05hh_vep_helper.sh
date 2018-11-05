#!/bin/bash
source ../../software/perlbrew_vep_environment.sh

set -euo pipefail

../../software/ensembl-vep/vep \
  --format vcf -i tmp/05_mgrb_split_variants.vcf.bgz \
  --tab -o tmp/05_mgrb_split_variants.vep.tab --no_stats --force_overwrite \
  --assembly GRCh37 --minimal --everything --pick --gencode_basic \
  --cache --offline --dir ../../software/.vep --fasta ../../software/.vep/homo_sapiens/90_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
  --fork 28
