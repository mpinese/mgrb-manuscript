import hail

hc = hail.HailContext(log = 'logs/05_rf_qualscore.log')

MGRB_LOCATION = '..'
GIAB_LOCATION = '../../GiaB_hail'

# Load the MGRB variants
vds_mgrb = hc.read('%s/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.vds' % MGRB_LOCATION)

# Create a consistent variant key in the MGRB data, to use to map the
# concordance results back to the full MGRB set.
vds_mgrb = vds_mgrb.annotate_variants_expr('va.variant_key = v.contig + ":" + v.start + ":" + v.ref + ":" + v.altAlleles.map(a => a.alt).mkString(",")')

# Extract the NA12878 variants only (ZAAAA).
import os.path
if os.path.exists('%s/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.ZAAAA.split.minrep.vds' % MGRB_LOCATION):
    vds_mgrb_na12878 = hc.read('%s/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.ZAAAA.split.minrep.vds' % MGRB_LOCATION)
else:
    vds_mgrb_na12878 = (vds_mgrb
        .filter_variants_table(hail.KeyTable.import_bed('%s/code/source_data/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed' % GIAB_LOCATION))
        .filter_samples_list(['ZAAAA'])
        .repartition(100)
        .split_multi()
        .min_rep()
    )
    vds_mgrb_na12878.write('%s/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.ZAAAA.split.minrep.vds' % MGRB_LOCATION)

print(vds_mgrb_na12878.summarize())


# Load the GiaB NA12878 variants
vds_giab_na12878 = (hc
    .read('%s/GiaB_HG001_3.3.2_highconf_nosomaticdel.vds' % GIAB_LOCATION)
    .split_multi()
    .min_rep()
    .rename_samples({'HG001': 'ZAAAA'})
)

summary, samples, variants = vds_mgrb_na12878.concordance(vds_giab_na12878)

# Classify MGRB ZAAAA variants as either correct or errors, contingent
# on a call having been made in MGRB.
# 
# Use the following confusion matrix:
#
#                 GiaB
#                   NoData  NoCall  HomRef     Het   HomVar
# MGRB    NoData        SZ      SZ      SZ      FN      FN
#         NoCall        TN      SZ      SZ      FN      FN
#         HomRef        TN      SZ      SZ      FN      FN
#         Het           FP      SZ      SZ      TP      DZ
#         HomVar        FP      SZ      SZ      DZ      TP
# 
# TP  True positive
# FP  False positive
# TN  True negative
# FN  False negative
# SZ  Structural zero
# DZ  Discordant zygosity
#
# Counts as of code writing:
# 
#                 GiaB                
#                   NoData  NoCall  HomRef     Het  HomVar
# MGRB    NoData         0       0       0    1106     818
#         NoCall     40304       0       0     210     191
#         HomRef  61645477       0       0    3147      83
#         Het       119500       0       0 2265390   26591
#         HomVar       954       0       0     252 1427515
#
# Rows and columns correspond to the summary list of lists, in row-major order.
#
# Construct an approximator to estimate Pr(Pos | U(Pos,Neg)), with sets
#   Pos = U(Het:Het,HomVar:HomVar), 
#   Neg = U(Het:NoData, HomVar:NoData, Het:HomVar, HomVar:Het)
# Where pairs are MGRB:GiaB
# 
# This is equivalent to Pr(Pos | U(Het:*, HomVar:*), or the probabilty
# of a Het or HomVar variant called in MGRB being correct.
#

variants_train = (variants
    .annotate('''
        pos = concordance[3][3] + concordance[4][4],
        neg = concordance[3][0] + concordance[4][0] + concordance[3][4] + concordance[4][3]''')
    .filter('pos != 0 || neg != 0')
    .select(['v', 'pos', 'neg'])
    .key_by('v')
)

# Annotate the MGRB NA12878 VDS with these training labels.
# This then enables the keying of the labels to the original
# full MGRB VDS variant IDs.

# Annotate the MGRB VDS with these training labels
vds_mgrb_na12878 = vds_mgrb_na12878.annotate_variants_table(variants_train, expr='''
    va.rf_var_train = {
        classVariantPos: table.pos,
        classVariantNeg: table.neg 
    }''')


vds_mgrb_na12878 = (vds_mgrb_na12878
    .annotate_variants_expr('''va.rf_var_train = {
        variantNumAlleles: v.nAlleles(),
        variantMixedType: v.altAlleles.map(a => a.category()).toSet().size() > 1,
        variantSpanningDel: v.altAlleles.exists(a => a.isStar()),
        variantInbreedingCoef: va.gatk.InbreedingCoeff,
        variantMQRankSum: va.gatk.MQRankSum,
        variantReadPosRankSum: va.gatk.ReadPosRankSum,
        variantSOR: va.gatk.SOR,
        variantVQSLOD: va.gatk.VQSLOD,
        variantQD: va.gatk.QD,
        variantHaplotypeScore: va.gatk.HaplotypeScore,
        variantBaseQRankSum: va.gatk.BaseQRankSum,
        variantClippingRankSum: va.gatk.ClippingRankSum,
        locusBadCoverage_dist: va.locus.badCoverage_dist,
        locusBadComplexity_dist: va.locus.badComplexity_dist,
        locusBadRepeat_dist: va.locus.badRepeat_dist,
        locusBadMappability_dist: va.locus.badMappability_dist,
        locusBadEncodeExcluded_dist: va.locus.badEncodeExcluded_dist,
        locusBadPAR: va.locus.badPAR,
        locusTier: va.locus.tier
    }''')
    .annotate_alleles_expr('''
        va.rf_var_train.classVariantPos = 0,
        va.rf_var_train.classVariantNeg = 0,
        va.rf_var_train.alleleType = v.altAllele.category(),
        va.rf_var_train.alleleMedianDP = gs.filter(g => g.isCalledNonRef()).map(g => g.fractionReadsRef()*g.dp).collect().median(),
        va.rf_var_train.alleleMedianGQ = gs.filter(g => g.isCalledNonRef()).map(g => g.gq).collect().median(),
        va.rf_var_train.alleleMedianDALT = gs.filter(g => g.isCalledNonRef()).map(g => g.dosage).collect().median(),
        va.rf_var_train.alleleMedianpAB = gs.filter(g => g.isHet()).map(g => g.pAB()).collect().median()
    ''')
)
vds_mgrb_na12878.summarize()


vds_mgrb = (vds_mgrb
    .annotate_variants_expr('''va.rf_var_train = {
        variantNumAlleles: v.nAlleles(),
        variantMixedType: v.altAlleles.map(a => a.category()).toSet().size() > 1,
        variantSpanningDel: v.altAlleles.exists(a => a.isStar()),
        variantInbreedingCoef: va.gatk.InbreedingCoeff,
        variantMQRankSum: va.gatk.MQRankSum,
        variantReadPosRankSum: va.gatk.ReadPosRankSum,
        variantSOR: va.gatk.SOR,
        variantVQSLOD: va.gatk.VQSLOD,
        variantQD: va.gatk.QD,
        variantHaplotypeScore: va.gatk.HaplotypeScore,
        variantBaseQRankSum: va.gatk.BaseQRankSum,
        variantClippingRankSum: va.gatk.ClippingRankSum,
        locusBadCoverage_dist: va.locus.badCoverage_dist,
        locusBadComplexity_dist: va.locus.badComplexity_dist,
        locusBadRepeat_dist: va.locus.badRepeat_dist,
        locusBadMappability_dist: va.locus.badMappability_dist,
        locusBadEncodeExcluded_dist: va.locus.badEncodeExcluded_dist,
        locusBadPAR: va.locus.badPAR,
        locusTier: va.locus.tier
    }''')
    .annotate_alleles_expr('''
        va.rf_var_train.classVariantPos = 0,
        va.rf_var_train.classVariantNeg = 0,
        va.rf_var_train.alleleType = v.altAllele.category(),
        va.rf_var_train.alleleMedianDP = gs.filter(g => g.isCalledNonRef()).map(g => g.fractionReadsRef()*g.dp).collect().median(),
        va.rf_var_train.alleleMedianGQ = gs.filter(g => g.isCalledNonRef()).map(g => g.gq).collect().median(),
        va.rf_var_train.alleleMedianDALT = gs.filter(g => g.isCalledNonRef()).map(g => g.dosage).collect().median(),
        va.rf_var_train.alleleMedianpAB = gs.filter(g => g.isHet()).map(g => g.pAB()).collect().median()
    ''')
)

#         NEGATIVE_TRAIN_SITE: Boolean,
#         POSITIVE_TRAIN_SITE: Boolean,

