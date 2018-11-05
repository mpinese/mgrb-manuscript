#!./bin/pyhail.sh
import hail

hc = hail.HailContext(log = 'log/01_gois.log', tmp_dir = 'tmp/hail')

#vds = hc.read('../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.vds')
vds = hc.read('temp.vds')

goi_kt = hail.KeyTable.import_bed('gois_refGene.bed')

(vds
#    .filter_variants_table(goi_kt)
#    .repartition(112)
#    .split_multi()
    .annotate_variants_expr('''va.info = {
        QUAL: va.variant.qual,
        tier: va.locus.tier,
        symbol: va.alleles[va.aIndex-1].annotation.vep.cross_references.symbol,
        consequences: va.alleles[va.aIndex-1].annotation.vep.consequences.consequences.mkString("|"),
        impact: va.alleles[va.aIndex-1].annotation.vep.consequences.impact_class,
        HGVSp: va.alleles[va.aIndex-1].annotation.vep.HGVS.HGVSp,
        MGRB_AC: va.alleles[va.aIndex-1].metrics.allele_counts.total,
        MGRB_AN: va.variant.metrics.chromatid_counts.total,
        MGRB_AF: va.alleles[va.aIndex-1].metrics.allele_frequencies.total,
        GnomAD_AF: va.alleles[va.aIndex-1].annotation.freqs.gnomad.AF}''')
    .filter_variants_expr('va.info.impact == "HIGH" || va.info.impact == "MODERATE"')
    .drop_samples()
    .export_variants('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusAnnot.WGStier12.gois.tsv', 'VARIANT = v, va.info.*')
#    .export_vcf('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusAnnot.WGStier12.gois.vcf.bgz')
)

