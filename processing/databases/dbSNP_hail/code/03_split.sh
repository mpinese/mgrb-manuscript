#!/bin/bash
set -euo pipefail

# Hail is failing to load the single-file dbSNP VCF.
# Try to help it along by splitting by chromosome.
parallel bcftools view -O z -o ../normed_data/All_20170710.normed.chrom{}.vcf.bgz ../normed_data/All_20170710.normed.vcf.bgz {} ::: $(seq 1 22) X Y MT
parallel tabix {} ::: ../normed_data/All_20170710.normed.chrom*.vcf.bgz
rm ../normed_data/All_20170710.normed.vcf.bgz
