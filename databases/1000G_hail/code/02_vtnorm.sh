#!/bin/bash

# In contrast to the MGRB data, in which min_rep only in Hail was used, we
# do perform VT norm on the 1000 genomes data, to bring them as much as is
# practical into alignment with the MGRB GATK output.

REFERENCE="../../../resources/hs37d5x/reference_genome/hs37d5x.fa"

mkdir -p ../phase3.normed.vcf

parallel ./bin/vt normalize -o ../phase3.normed.vcf/{/} -r "${REFERENCE}" {} ::: ../phase3.vcf/*.vcf.gz 2>&1 | tee -a log/02_vtnorm.log
parallel mv {} {.}.bgz ::: ../phase3.normed.vcf/*.gz
parallel tabix {} ::: ../phase3.normed.vcf/*.bgz

# Hail doesn't handle the 1000G coding of Y genotypes (in which female data are missing) --
# rename the chrY files to effectively suppress loading by hail.
mv ../phase3.normed.vcf/ALL.chrY.phase3_integrated_v1b.20130502.genotypes.vcf.bgz ../phase3.normed.vcf/ALL.chrY.phase3_integrated_v1b.20130502.genotypes.vcf.gz
mv ../phase3.normed.vcf/ALL.chrY.phase3_integrated_v1b.20130502.genotypes.vcf.bgz.tbi ../phase3.normed.vcf/ALL.chrY.phase3_integrated_v1b.20130502.genotypes.vcf.gz.tbi


