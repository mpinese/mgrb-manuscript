#!/bin/bash
set -euo pipefail

REFERENCE="../../../resources/hs37d5x/reference_genome/hs37d5x.fa"

mkdir -p ../normed_data log

./bin/vt normalize -o ../normed_data/All_20170710.normed.vcf.gz -r "${REFERENCE}" ../source_data/All_20170710.vcf.gz 2>&1 | tee -a log/02_vtnorm.log

mv ../normed_data/All_20170710.normed.vcf.gz ../normed_data/All_20170710.normed.vcf.bgz
tabix ../normed_data/All_20170710.normed.vcf.bgz

