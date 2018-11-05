#!/bin/bash
set -euo pipefail
[ ! -e ../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusAnnot.WGStier12.ACMG.vcf.gz ] && ln -s MGRB.phase2.SNPtier12.match.vqsr.minrep.locusAnnot.WGStier12.ACMG.vcf.bgz ../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusAnnot.WGStier12.ACMG.vcf.gz
[ ! -e ../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusAnnot.WGStier12.ACMG.vcf.gz.tbi ] && tabix ../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusAnnot.WGStier12.ACMG.vcf.gz

../../../software/ensembl-vep/vep \
  --fork 7 -v \
  -i ../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusAnnot.WGStier12.ACMG.vcf.gz --format vcf \
  -o ../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusAnnot.WGStier12.ACMG.VEP89_noCommon_gencodebasic.vcf.gz --vcf --compress_output bgzip \
  --everything --hgvsg --pick --no_escape \
  --filter_common --no_intergenic --gencode_basic \
  --offline --dir_cache ../../software/.vep/ --fasta ../../../resources/hs37d5x/reference_genome/hs37d5x.fa

gzip -dc ../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusAnnot.WGStier12.ACMG.VEP89_noCommon_gencodebasic.vcf.gz | \
  grep -E '(^#)|MODERATE|HIGH' | \
  bgzip > ../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusAnnot.WGStier12.ACMG.VEP89_noCommon_gencodebasic.moderatehigh.vcf.gz

gzip -dc ../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusAnnot.WGStier12.ACMG.VEP89_noCommon_gencodebasic.moderatehigh.vcf.gz | \
  python3 bin/vcf2long_simple.py > ../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusAnnot.WGStier12.ACMG.VEP89_noCommon_gencodebasic.moderatehigh.tsv
