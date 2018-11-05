#!./bin/pyhail.sh
import hail
import sys
import pprint

hc = hail.HailContext(tmp_dir = 'tmp/hail')
vds = hc.read('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.vds')

vds.export_samples('99_samples.tsv', 'SampleID = s, tier1 = sa.tier1, tier2 = sa.tier2, pheno = sa.pheno.*')

