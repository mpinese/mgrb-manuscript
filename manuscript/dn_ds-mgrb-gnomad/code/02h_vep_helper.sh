#!/bin/bash
source ../../../software/perlbrew_vep_environment.sh

set -euo pipefail

#  --coding_only --gencode_basic --chr 1-22 \

../../../software/ensembl-vep/vep \
  --coding_only --gencode_basic \
  --format vcf -i ../02_gnomad.genomes.r2.0.1.sites.autosomes.split.all.freqs.minimal.vcf.gz \
  --tab -o ../02_gnomad.genomes.r2.0.1.sites.autosomes.split.all.freqs.minimal.vep.tab --no_stats --force_overwrite \
  --assembly GRCh37 --minimal --everything --pick --gencode_basic \
  --cache --offline --dir ../../../software/.vep --fasta ../../../software/.vep/homo_sapiens/90_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
  --fork 28

../../../software/ensembl-vep/vep \
  --coding_only --gencode_basic \
  --format vcf -i ../02_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.all.freqs.minimal.vcf.gz \
  --tab -o ../02_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.all.freqs.minimal.vep.tab --no_stats --force_overwrite \
  --assembly GRCh37 --minimal --everything --pick --gencode_basic \
  --cache --offline --dir ../../../software/.vep --fasta ../../../software/.vep/homo_sapiens/90_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
  --fork 28

