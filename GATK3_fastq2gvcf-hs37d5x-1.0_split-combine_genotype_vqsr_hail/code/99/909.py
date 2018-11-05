#!./bin/pyhail.sh
import hail

hc = hail.HailContext(tmp_dir = 'tmp/hail')

vds = hc.read('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.vds')

vds.drop_samples().annotate_variants_expr('va = select(va, locus)').write('temp_locus.vds')
