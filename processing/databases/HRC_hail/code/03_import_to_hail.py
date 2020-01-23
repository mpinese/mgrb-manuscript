#!./bin/pyhail.sh
import hail

hc = hail.HailContext(log = 'log/03_import_to_hail.log', tmp_dir = 'tmp/hail')

vds = hc.import_vcf('../normed_data/HRC.r1-1.GRCh37.wgs.mac5.sites.normed.vcf.bgz')

vds = (vds
    .annotate_variants_expr('''va.info = {
        AC: va.info.AC[0],
        AN: va.info.AN,
        AF: va.info.AF[0],
        AC_EXCLUDING_1000G: va.info.AC_EXCLUDING_1000G[0],
        AN_EXCLUDING_1000G: va.info.AN_EXCLUDING_1000G,
        AF_EXCLUDING_1000G: va.info.AF_EXCLUDING_1000G[0]}''')
    .min_rep()
    .repartition(560).write('../HRC.minrep.vds')
)
