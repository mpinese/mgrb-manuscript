#!./bin/pyhail.sh
import hail

hc = hail.HailContext(log = 'log/03_import_to_hail.log', tmp_dir = 'tmp/hail')

vds = hc.import_vcf('../clinvar.vcf.bgz').repartition(28).min_rep().write('../clinvar.vds')
