#!bin/pyhail.sh
import hail
#import pprint

hc = hail.HailContext(log = 'log/04_extract_na12878.log', tmp_dir = 'tmp/hail')

vds = hc.read('../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.vds')

#pprint.pprint(vds.variant_schema)

vds_control = (vds
    .filter_samples_expr('s == "ZAAAA"')
    .filter_variants_expr('gs.filter(g => g.isCalledNonRef()).count() > 0')
    .annotate_variants_expr('va = {info: {tier: va.locus.tier } }')
    .set_va_attributes('va.info.tier', {'Description': 'Locus confidence tier', 'Number': '1'})
)

vds_control.export_vcf('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.NA12878.vcf.bgz')
