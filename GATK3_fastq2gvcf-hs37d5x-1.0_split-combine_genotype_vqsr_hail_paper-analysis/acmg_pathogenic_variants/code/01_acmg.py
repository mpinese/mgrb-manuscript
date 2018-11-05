#!./bin/pyhail.sh
import hail

hc = hail.HailContext(log = 'log/01_acmg.log', tmp_dir = 'tmp/hail')

vds = hc.read('../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.vds')

acmg_kt = hail.KeyTable.import_bed('ACMG_refGene.bed')

vds = (vds
    .filter_variants_table(acmg_kt)
    .repartition(112)
    .filter_variants_expr('gs.filter(g => g.isCalledNonRef).count() > 0')
    .annotate_variants_expr('''va.info = {
        callRate: gs.fraction(g => g.isCalled),
        isSingleton: gs.filter(g => g.isCalledNonRef).count() == 1,
        NS: gs.filter(g => g.isCalled).count()}''')
    .annotate_alleles_expr('va.info.AC = gs.filter(g => g.isHet).count() + gs.filter(g => g.isHomVar).count()*2')
    .annotate_variants_expr('va.info.AF = va.info.AC / (va.info.NS*2)')
    .filter_alleles('va.info.AC[aIndex - 1] == 0', annotation='''
        va.info.AC = aIndices[1:].map(i => va.info.AC[i - 1]),
        va.info.AF = aIndices[1:].map(i => va.info.AF[i - 1])''',
        keep=False, subset=True, keep_star=False)
    .annotate_variants_expr('''va.info = {
        tier: va.locus.tier, 
        badCoverage_dist: va.locus.badCoverage_dist,
        badComplexity_dist: va.locus.badComplexity_dist, 
        badRepeat_dist: va.locus.badRepeat_dist, 
        badMappability_dist: va.locus.badMappability_dist,
        badEncodeExcluded_dist: va.locus.badEncodeExcluded_dist,
        goodGiaBHighConf_dist: va.locus.goodGiaBHighConf_dist,
        badPAR: va.locus.badPAR}''')
    .write('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusAnnot.WGStier12.ACMG.vds')
)

#vds = hc.read('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusAnnot.WGStier12.ACMG.vds')

#        if (!(v.contig ~ "^([0-9]+|X|Y)$")) 3

vds = (vds
    .annotate_variants_expr('''va.locus.tier =
        (if (va.locus.badCoverage_dist <= 5 || va.locus.badComplexity_dist <= 5 || va.locus.badMappability_dist <= 5 || va.locus.badEncodeExcluded_dist <= 5 || va.locus.badPAR) 3
        else (if (va.locus.goodGiaBHighConf_dist > 0) 2
        else (if (v.contig ~ "^(X|Y)$") 2
        else 1 )))''')
    .annotate_variants_expr('va.info.tier = va.locus.tier')
)

vds.export_vcf('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusAnnot.WGStier12.ACMG.vcf.bgz')
