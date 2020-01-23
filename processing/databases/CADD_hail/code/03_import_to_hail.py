#!./bin/pyhail.sh
import os
import os.path
import shutil
import hail

hc = hail.HailContext(log = 'log/03_import_to_hail.log', tmp_dir = 'tmp/hail')

hc.import_vcf('../converted_data/*.vcf.bgz', min_partitions=20000).repartition(10000).min_rep().write('../CADD_1.3.vds')
