#!./bin/pyhail.sh
import hail
import sys
import pprint

hc = hail.HailContext(tmp_dir = 'tmp/hail')
vds = hc.read('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.vds')

# Exclude control samples, alleles without sample support, and star alleles from summary
vds = vds.filter_samples_expr('"^Z" ~ s', keep=False).filter_alleles('va.alleles[aIndex - 1].metrics.allele_frequencies.total == 0 || v.altAlleles[aIndex - 1].isStar()', keep=False)

pprint.pprint(vds.summarize())

print('Singletons: {}'.format(vds.filter_variants_expr('gs.filter(g => g.isCalledNonRef()).count() == 1').count_variants()))
