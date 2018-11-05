#!./bin/pyhail.sh
import hail

hc = hail.HailContext(log = 'log/01h_export_freqs_gnomadwgs_all.log', tmp_dir = 'tmp/hail')

vds = hc.read('../gnomad.genomes.r2.0.1.sites.autosomes.vds')

(vds
    .filter_variants_expr('v.isAutosomal && va.filters.isEmpty() && va.pass')
    .split_multi()
    .filter_variants_expr('v.altAllele().isSNP()')
    .annotate_variants_expr('''va = {info: {
        AF_NFE: va.info.AF_NFE[va.aIndex-1],
        AC_NFE: va.info.AC_NFE[va.aIndex-1],
        AN_NFE: va.info.AN_NFE,
        nHomVar_NFE: va.info.Hom_NFE[va.aIndex-1]}}''')
    .set_va_attributes('va.info.AF_NFE', {'Description': 'Calculated alternate allele frequency, GnomAD NFE', 'Number': '1'})
    .set_va_attributes('va.info.AC_NFE', {'Description': 'Count of alternate alleles, GnomAD NFE', 'Number': '1'})
    .set_va_attributes('va.info.AN_NFE', {'Description': 'Count of alleles, GnomAD NFE', 'Number': '1'})
    .set_va_attributes('va.info.nHomVar_NFE', {'Description': 'Number of homozygous alternate samples, GnomAD NFE', 'Number': '1'})
    .export_vcf('../01_gnomad.genomes.r2.0.1.sites.autosomes.split.all.freqs.vcf.bgz')
)

#vds.export_variants('../01_gnomad.genomes.r2.0.1.sites.autosomes.split.goi.all.freqs.tsv', 'variant=v, va.info.*')

