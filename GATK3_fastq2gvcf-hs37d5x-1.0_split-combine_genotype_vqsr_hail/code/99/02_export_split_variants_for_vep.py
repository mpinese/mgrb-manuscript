#!../../software/pyhail.sh
import hail
from hail.expr import TString, TBoolean

hc = hail.HailContext(log = 'log/02_export_split_variants_for_vep.log', tmp_dir = 'tmp/hail')

vds = hc.read('../MGRB.phase2.tier12.match.vqsr.minrep.vds')


vds = vds.split_multi().drop_samples().repartition(28, shuffle=False).annotate_variants_expr(['va.rsid = str(v)', 'va.info = {}'])
vds.export_vcf('../MGRB.phase2.tier12.match.vqsr.minrep.split.variants.vcf.bgz', parallel=True)
