#!/bin/bash
REFERENCE="/g/data3/wq2/resources/hs37d5x/reference_genome/hs37d5x.fa"

ls -1 ../*.vcf.gz > input.list
java -cp bin/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants -R "${REFERENCE}" -V input.list -out ../MGRB.phase2qc.hc.vcf.gz -assumeSorted

./bin/vt normalize -r "${REFERENCE}" -o ../MGRB.phase2qc.hc.norm.vcf.gz ../MGRB.phase2qc.hc.vcf.gz
tabix ../MGRB.phase2qc.hc.norm.vcf.gz

ln -s MGRB.phase2qc.hc.norm.vcf.gz ../MGRB.phase2qc.hc.norm.vcf.bgz
ln -s MGRB.phase2qc.hc.norm.vcf.gz.tbi ../MGRB.phase2qc.hc.norm.vcf.bgz.tbi

