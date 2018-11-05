#!./bin/pyhail.sh
import hail

hc = hail.HailContext(log = 'log/02_split.log', tmp_dir = 'tmp/hail')

vds_gnomad_x = hc.read('../gnomad.genomes.r2.0.1.sites.X.vds')
vds_gnomad_autosomes = hc.read('../gnomad.genomes.r2.0.1.sites.autosomes.vds')

(vds_gnomad_x
    .annotate_variants_expr('''va.info = select(va.info,
        AF_AFR, AF_AMR, AF_ASJ, AF_EAS, AF_FIN, AF_NFE, AF_OTH, AF,
        AC_AFR, AC_AMR, AC_ASJ, AC_EAS, AC_FIN, AC_NFE, AC_OTH, AC,
        AN_AFR, AN_AMR, AN_ASJ, AN_EAS, AN_FIN, AN_NFE, AN_OTH, AN)''')
    .annotate_variants_expr('va = select(va, info)')
    .write('../gnomad.genomes.r2.0.1.sites.X.simplified.vds')
)

(vds_gnomad_autosomes
    .annotate_variants_expr('''va.info = select(va.info,
        AF_AFR, AF_AMR, AF_ASJ, AF_EAS, AF_FIN, AF_NFE, AF_OTH, AF,
        AC_AFR, AC_AMR, AC_ASJ, AC_EAS, AC_FIN, AC_NFE, AC_OTH, AC,
        AN_AFR, AN_AMR, AN_ASJ, AN_EAS, AN_FIN, AN_NFE, AN_OTH, AN)''')
    .annotate_variants_expr('va = select(va, info)')
    .write('../gnomad.genomes.r2.0.1.sites.autosomes.simplified.vds')
)
