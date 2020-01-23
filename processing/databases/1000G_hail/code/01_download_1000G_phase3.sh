#!/bin/bash
set -euo pipefail

mkdir -p ../phase3.vcf
cd ../phase3.vcf
parallel --gnu wget -nc ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ::: $(seq 1 22)
parallel --gnu wget -nc ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi ::: $(seq 1 22)
wget -nc ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chrY.phase3_integrated_v1b.20130502.genotypes.vcf.gz&
wget -nc ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chrY.phase3_integrated_v1b.20130502.genotypes.vcf.gz.tbi&
wget -nc ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz&
wget -nc ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz.tbi&
wget -nc ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/README_populations.md&
curl ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel | sed -E 's/\t+$//' > integrated_call_samples_v3.20130502.ALL.panel

wait

cd -

