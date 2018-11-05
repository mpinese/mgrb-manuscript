#!./bin/pyhail.sh
import hail

hc = hail.HailContext(log = 'log/01_export_nfe_tab.log', tmp_dir = 'tmp/hail')

vds_gnomad_split = hc.read('../../GnomAD_hail/gnomad.genomes.r2.0.1.sites.combined.split.minrep.vds')
vds_dbsnp = hc.read('../../dbSNP_hail/dbSNP_150.split.minrep.vds')

(vds_gnomad_split
    .filter_variants_expr('v.isAutosomal()')
    .annotate_variants_expr('va.info = {AF_NFE: va.gnomad.AF_NFE, AC_NFE: va.gnomad.AC_NFE, AN_NFE: va.gnomad.AN_NFE}')
    .filter_variants_expr('va.info.AC_NFE > 0 && va.info.AN_NFE >= 5000')
    .annotate_variants_vds(vds_dbsnp, expr='va.rsid = vds.info.rsid')
    .export_variants('../01_gnomad.genomes.r2.0.1.sites.combined.split.minrep.NFE.dbSNP.autosomes.tsv', 'v.contig, v.start, v.ref, v.alt(), va.rsid, va.info.AC_NFE, va.info.AN_NFE')
)

#    .set_va_attributes('va.info.AF_NFE', {'Description': 'NFE alternate allele frequency', 'Number': '1'})
#    .set_va_attributes('va.info.AC_NFE', {'Description': 'NFE alternate allele count', 'Number': '1'})
#    .set_va_attributes('va.info.AN_NFE', {'Description': 'NFE total allele count', 'Number': '1'})
#    .export_vcf('../01_gnomad.genomes.r2.0.1.sites.combined.split.minrep.NFE.dbSNP.autosomes.vcf.bgz')
