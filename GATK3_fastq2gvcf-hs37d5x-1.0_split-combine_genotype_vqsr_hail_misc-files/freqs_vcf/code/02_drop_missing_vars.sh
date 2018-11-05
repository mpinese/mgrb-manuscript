#!/bin/sh

gzip -dc ./tmp/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.freqs.vcf.bgz | grep -vF 'AC=0;' | bgzip > ../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.freqs.vcf.bgz
