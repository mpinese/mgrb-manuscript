#!./bin/pyhail.sh
import hail
import sys
import pprint

hc = hail.HailContext(tmp_dir = 'tmp/hail')
vds = hc.read('../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.vds')

vds_notcontrol = vds.filter_samples_expr('"^[^Z]" ~ s').filter_variants_expr('gs.filter(g => g.isCalledNonRef()).count() > 0')
vds_aspree = vds.filter_samples_expr('"^A" ~ s').filter_variants_expr('gs.filter(g => g.isCalledNonRef()).count() > 0')
vds_45nup = vds.filter_samples_expr('"^B" ~ s').filter_variants_expr('gs.filter(g => g.isCalledNonRef()).count() > 0')
vds_control = vds.filter_samples_expr('"^Z" ~ s').filter_variants_expr('gs.filter(g => g.isCalledNonRef()).count() > 0')

print('Total:')
pprint.pprint(vds.summarize())
print()
print('Not control:')
pprint.pprint(vds_notcontrol.summarize())
print()
print('ASPREE:')
pprint.pprint(vds_aspree.summarize())
print()
print('45nUp:')
pprint.pprint(vds_45nup.summarize())
print()
print('Control:')
pprint.pprint(vds_control.summarize())
print()
