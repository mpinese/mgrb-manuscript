#!/bin/bash
set -euo pipefail

curl https://www.ebi.ac.uk/gwas/api/search/downloads/alternative -o gwas_catalogue.tsv

wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz.tbi -O dbsnp_All_20180423.vcf.gz.tbi
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz -O dbsnp_All_20180423.vcf.gz
