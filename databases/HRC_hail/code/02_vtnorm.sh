#!/bin/bash
set -euo pipefail

# In contrast to the MGRB data, in which min_rep only in Hail was used, we
# do perform VT norm on the 1000 genomes data, to bring them as much as is
# practical into alignment with the MGRB GATK output.

REFERENCE="../../../resources/hs37d5x/reference_genome/hs37d5x.fa"

mkdir -p ../normed_data

./bin/vt normalize -o ../normed_data/HRC.r1-1.GRCh37.wgs.mac5.sites.normed.vcf.gz -r "${REFERENCE}" ../source_data/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz | tee -a log/02_vtnorm.log

mv ../normed_data/HRC.r1-1.GRCh37.wgs.mac5.sites.normed.vcf.gz ../normed_data/HRC.r1-1.GRCh37.wgs.mac5.sites.normed.vcf.bgz
tabix ../normed_data/HRC.r1-1.GRCh37.wgs.mac5.sites.normed.vcf.bgz

