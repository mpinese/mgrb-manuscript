#!/bin/bash

gzip -dc ../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.NA12878.vcf.gz | \
    python ./05_drop_missing_alleles.py | bgzip > ../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.NA12878.only.vcf.gz

gzip -dc ../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.NA12878.only.vcf.gz | \
    grep -vF 'tier=3' | bgzip > ../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.NA12878.only.tier12.vcf.gz
