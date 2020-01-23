#!./bin/pyhail.sh
import hail

hc = hail.HailContext(tmp_dir = 'tmp/hail')

hc.import_vcf('test_vars.vcf.bgz').write('test_vars.vds')
hc.import_vcf('test_annots.vcf.bgz').write('test_annots.vds')


