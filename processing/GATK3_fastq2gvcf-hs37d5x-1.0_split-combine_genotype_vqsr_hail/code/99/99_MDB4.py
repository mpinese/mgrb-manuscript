#!./bin/pyhail.sh
import hail
from hail.representation import Interval

hc = hail.HailContext(tmp_dir = 'tmp/hail')

#vds = hc.read('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.vds/')
#
#vds = (vds
#    .filter_intervals(Interval.parse('3:129140000-129160000'))
#    .filter_variants_expr('va.alleles.exists(e => e.annotation.vep.consequences.impact_class == "HIGH")')
#)
#
#vds.write('test.vds')

vds = hc.read('test.vds')
vds.export_genotypes('test.tsv', 'SAMPLE=s, VARIANT=v, CSQ=va.alleles.map(a => a.annotation.vep.consequences.consequences.mkString(",")).mkString(";"), GENOTYPE=g')
