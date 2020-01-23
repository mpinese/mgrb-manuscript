#!./bin/pyhail.sh
import hail

hc = hail.HailContext(log = 'log/03_import_to_hail.log', tmp_dir = 'tmp/hail')

(hc
    .import_vcf('../converted_data/Eigen_hg19_coding_annot_04092016.norm.vcf.bgz')
    .min_rep()
    .deduplicate()
    .write('../Eigen_coding_04092016.vds')
)


