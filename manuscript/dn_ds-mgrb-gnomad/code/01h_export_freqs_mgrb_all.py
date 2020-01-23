#!./bin/pyhail.sh
import hail

hc = hail.HailContext(log = 'log/01h_export_goi_freqs_mgrb_all.log', tmp_dir = 'tmp/hail')

vds = hc.read('../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.vds')

vds = (vds
    .filter_variants_expr('v.isAutosomal && va.variant.filters.isEmpty() && va.locus.tier != 3')
    .split_multi()
    .filter_variants_expr('v.altAllele().isSNP()')
    .variant_qc()
    .drop_samples()
    .annotate_variants_expr('va = {info: {AF_MGRB: va.qc.AF, AC_MGRB: va.qc.AC, AN_MGRB: va.qc.nCalled*2, nHomVar_MGRB: va.qc.nHomVar}}')
    .set_va_attributes('va.info.AF_MGRB', {'Description': 'Calculated alternate allele frequency, MGRB', 'Number': '1'})
    .set_va_attributes('va.info.AC_MGRB', {'Description': 'Count of alternate alleles, MGRB', 'Number': '1'})
    .set_va_attributes('va.info.AN_MGRB', {'Description': 'Count of alleles, MGRB', 'Number': '1'})
    .set_va_attributes('va.info.nHomVar_MGRB', {'Description': 'Number of homozygous alternate samples, MGRB', 'Number': '1'})
)

vds.export_vcf('../01_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.all.freqs.vcf.bgz')

#vds.export_variants('../01_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.goi.all.freqs.tsv', 'variant=v, va.info.*')
