#!/bin/bash
set -euo pipefail

mkdir -p ../source_data
cd ../source_data
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/All_20170710.vcf.gz
tabix All_20170710.vcf.gz
date > dbsnp.download_date
cd -
