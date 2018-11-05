#!/bin/bash
set -euo pipefail

rm -f ../*.vqsr_snp.vcf.gz
rm -f ../*.vqsr_snp.vcf.gz.tbi

for f in ../*.vcf.gz.vqsr_snp.vcf.gz.vqsr_indel.vcf.gz; do
  mv ${f} ${f%.vcf.gz.vqsr_snp.vcf.gz.vqsr_indel.vcf.gz}.vqsr.vcf.bgz
  mv ${f}.tbi ${f%.vcf.gz.vqsr_snp.vcf.gz.vqsr_indel.vcf.gz}.vqsr.vcf.bgz.tbi
done

