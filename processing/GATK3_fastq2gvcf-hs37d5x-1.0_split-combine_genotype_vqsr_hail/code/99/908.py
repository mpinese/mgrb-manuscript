#!./bin/pyhail.sh
import hail

hc = hail.HailContext(tmp_dir = 'tmp/hail')

vds = hc.read('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.nolocus.vds')
loc_vds = hc.read('temp_locus.vds')

vds = vds.annotate_variants_vds(loc_vds, expr='va.locus = vds.locus')

vds.write('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.vds')
