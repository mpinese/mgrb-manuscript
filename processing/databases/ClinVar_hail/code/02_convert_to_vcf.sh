#!/bin/bash
set -euo pipefail

mkdir -p tmp

python2 clinvar2vcf.py ../source_data/variant_summary.txt.gz ../../../resources/hs37d5x/reference_genome/hs37d5x.fa > tmp/clinvar.unsorted.vcf

head -n 1000 tmp/clinvar.unsorted.vcf | grep '^#' > tmp/header.txt

grep -v '^#' tmp/clinvar.unsorted.vcf | LC_ALL=C sort -k1,1 -k2,2n > tmp/clinvar_sorted.noheader

cat tmp/header.txt tmp/clinvar_sorted.noheader | bgzip > tmp/clinvar.vcf.gz
tabix tmp/clinvar.vcf.gz

./bin/vt normalize -o ../clinvar.vcf.gz -r ../../../resources/hs37d5x/reference_genome/hs37d5x.fa tmp/clinvar.vcf.gz
mv ../clinvar.vcf.gz ../clinvar.vcf.bgz
tabix ../clinvar.vcf.bgz
