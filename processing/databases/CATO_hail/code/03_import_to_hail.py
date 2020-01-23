#!./bin/pyhail.sh
import hail

hc = hail.HailContext(log = 'log/03_import_to_hail.log', tmp_dir = 'tmp/hail')

(hc
    .import_vcf('../converted_data/dbSNP142.CATO.V1.1.norm.vcf.bgz')
    .min_rep()
    .deduplicate()
    .write('../CATO_1.1.vds')
)
