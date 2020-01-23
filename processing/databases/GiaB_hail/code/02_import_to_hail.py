#!./bin/pyhail.sh
import hail

hc = hail.HailContext(log = 'log/02_import_to_hail.log', tmp_dir = 'tmp/hail')

vds = (hc
    .import_vcf('source_data/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.minimal.vcf.bgz')
    .filter_variants_table(hail.KeyTable.import_bed('source_data/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed'))
    .repartition(100)
)

vds.write('../GiaB_HG001_3.3.2_highconf_nosomaticdel.vds')
