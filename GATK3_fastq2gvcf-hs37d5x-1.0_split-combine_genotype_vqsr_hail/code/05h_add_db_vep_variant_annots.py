#!./bin/pyhail.sh
import hail
import os
import pprint
from hail.expr import TVariant



hc = hail.HailContext(log = 'log/05h_add_db_vep_variant_annots.log', tmp_dir = 'tmp/hail')

vds_mgrb = hc.read('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.vds').repartition(10000)

# Generate a temporary frequency-only version of the MGRB.  This is in part to
# act as a QC for the autosomal MGRB frequencies already calculated, as well as
# to provide sex naive frequency estimates for better comparison to other cohorts
# which don't do sex-aware AF calculation.
vds_mgrb.split_multi().annotate_variants_expr('va = {}').variant_qc().drop_samples().write('tmp/05_mgrb_simpleafs.vds')

# Annotate MGRB with the 1000G, HRC, and GnomAD variant frequencies.  These frequency
# databases are all split, but MGRB is multi-allelic.  As a consequence, some gymnastics 
# are needed to perform this annotation correctly.

# What we do is create a variant-only multi-allelic MGRB vds, then progressively add the
# split annotations to this vds.  Finally, this 'annotation' MGRB VDS is used to directly
# annotate the final MGRB files.

# Annotation functions by lfrancoli 
# (http://discuss.hail.is/t/annotating-a-non-split-vds-using-a-split-vds/123)
# Example usage:
# annotate_non_split_from_split(hc, non_split_vds_path='gs://mybucket/my.vds',
#                               split_vds=hc.read("gs://mybucket/split_vds"),
#                               annotations=['va.annotation1','va.foo.annotation2'])
def getAnnType(annotation, schema):
    ann_path = annotation.split(".")[1:]
    ann_type = schema
    for p in ann_path:
        ann_type = [x for x in ann_type.fields if x.name == p][0].typ
    return ann_type


def annotate_non_split_from_split(hc, non_split_vds_path, split_vds, annotations):
    ann_types = list(map(lambda x: str(getAnnType(x, split_vds.variant_schema)), annotations))

    variant_annotated_vds = (hc
        .read(non_split_vds_path, drop_samples=True)
        .annotate_variants_expr('va.variant = v, va.nAltAlleles = v.nAltAlleles()')
        .split_multi()
    )

    ann_agg_codes = ["`%s` = index(va.map(x => {val: %s, aIndex: va.aIndex}).collect(), aIndex)" % (a, a) for a in annotations]

    agg = (split_vds
        .annotate_variants_vds(variant_annotated_vds, 'va.variant = vds.variant, va.aIndex = vds.aIndex, va.nAltAlleles = vds.nAltAlleles')
        .filter_variants_expr('isDefined(va.variant)')
        .variants_table().repartition(10000)
        .aggregate_by_key('variant = va.variant', ann_agg_codes + ['nAltAlleles = va.map(x => x.nAltAlleles).take(1)[0]'])
    )

    ann_codes = ['%s = let x = table.`%s` in' \
                 ' range(table.nAltAlleles).map(i => if(x.contains(i+1)) x[i+1].val else NA: %s)' % (a, a, b)
                 for (a, b) in zip(annotations, ann_types)]

    return (hc
        .read(non_split_vds_path)
        .annotate_variants_table(agg, expr=",".join(ann_codes))
    )


# Create the MGRB annotation-only VDS starting point
#######################################################################################################################################################################
vds_mgrb.drop_samples().repartition(10000).write('tmp/05_mgrb_varsonly.vds')

# Load the source annotation databases

vds_1000g_freqs = hc.read('../../databases/1000G_hail/1000G.split.minrep.freqs.vds').repartition(10000)
vds_hrc_freqs = hc.read('../../databases/HRC_hail/HRC.minrep.vds').repartition(10000)
vds_gnomad_freqs = hc.read('../../databases/GnomAD_hail/gnomad.genomes.r2.0.1.sites.combined.split.minrep.vds').repartition(10000)
vds_mgrb_split_freqs = hc.read('tmp/05_mgrb_simpleafs.vds').repartition(10000)
vds_dbsnp = hc.read('../../databases/dbSNP_hail/dbSNP_150.split.minrep.vds').repartition(10000)
vds_clinvar = hc.read('../../databases/ClinVar_hail/clinvar.vds').repartition(10000)
vds_cato = hc.read('../../databases/CATO_hail/CATO_1.1.vds').repartition(10000)
vds_eigen = hc.read('../../databases/Eigen_hail/Eigen_coding_04092016.vds').repartition(10000)
vds_vep = hc.read('tmp/05_mgrb_split_variants.vep.vds')
#vds_cadd = hc.read('../../databases/CADD_hail/CADD_1.3.vds')

# Relabel the annotations in the source databases to avoid collisions
vds_1000g_freqs = vds_1000g_freqs.annotate_variants_expr('''
    va.freqs = {tgp: {
        AC: va.info.AC,
        AF: va.info.AF,
        NS: va.info.NS,
        AN: va.info.AN,
        AF_EAS: va.info.EAS_AF,
        AF_EUR: va.info.EUR_AF,
        AF_AFR: va.info.AFR_AF,
        AF_AMR: va.info.AMR_AF,
        AF_SAS: va.info.SAS_AF,
        AC_Hail: va.qc.AC,
        AF_Hail: va.qc.AF,
        NS_Hail: va.qc.nCalled
    }}''')

vds_hrc_freqs = vds_hrc_freqs.annotate_variants_expr('''
    va.freqs = {hrc: {
        AC: va.info.AC,
        AN: va.info.AN,
        AF: va.info.AF,
        AC_EXCLUDING_1000G: va.info.AC_EXCLUDING_1000G,
        AN_EXCLUDING_1000G: va.info.AN_EXCLUDING_1000G,
        AF_EXCLUDING_1000G: va.info.AF_EXCLUDING_1000G
    }}''')

vds_gnomad_freqs = vds_gnomad_freqs.annotate_variants_expr('''
    va.freqs = {gnomad: {
        AC: va.gnomad.AC,
        AF: va.gnomad.AF,
        AN: va.gnomad.AN,
        AC_AFR: va.gnomad.AC_AFR,
        AC_AMR: va.gnomad.AC_AMR,
        AC_ASJ: va.gnomad.AC_ASJ,
        AC_EAS: va.gnomad.AC_EAS,
        AC_FIN: va.gnomad.AC_FIN,
        AC_NFE: va.gnomad.AC_NFE,
        AC_OTH: va.gnomad.AC_OTH,
        AN_AFR: va.gnomad.AN_AFR,
        AN_AMR: va.gnomad.AN_AMR,
        AN_ASJ: va.gnomad.AN_ASJ,
        AN_EAS: va.gnomad.AN_EAS,
        AN_FIN: va.gnomad.AN_FIN,
        AN_NFE: va.gnomad.AN_NFE,
        AN_OTH: va.gnomad.AN_OTH,
        AF_AFR: va.gnomad.AF_AFR,
        AF_AMR: va.gnomad.AF_AMR,
        AF_ASJ: va.gnomad.AF_ASJ,
        AF_EAS: va.gnomad.AF_EAS,
        AF_FIN: va.gnomad.AF_FIN,
        AF_NFE: va.gnomad.AF_NFE,
        AF_OTH: va.gnomad.AF_OTH
    }}''')

vds_mgrb_split_freqs = vds_mgrb_split_freqs.annotate_variants_expr('''
    va.freqs = {mgrb: {
        AC_Hail: va.qc.AC,
        AF_Hail: va.qc.AF,
        NS_Hail: va.qc.nCalled
    }}''')

vds_clinvar = vds_clinvar.annotate_variants_expr('''
    va.clinvar = {
        AlleleID: va.info.AlleleID,
        ClinicalSignificance: va.info.ClinicalSignificance,
        ClinSigSimple: va.info.ClinSigSimple,
        rsid: va.info.rsid,
        ReviewStatus: va.info.ReviewStatus,
        NumberSubmitters: va.info.NumberSubmitters
    }''')

vds_cato = vds_cato.annotate_variants_expr('''
    va.predictions = {CATO: {
         score: va.info.score,
         motif: va.info.motif,
         strand: va.info.strand,
         motifpos: va.info.motifpos,
         celltypes: va.info.celltypes
    }}''')

vds_eigen = vds_eigen.annotate_variants_expr('''
    va.predictions = {Eigen: {
        EigenRaw: va.info.EigenRaw,
        EigenPhred: va.info.EigenPhred,
        EigenPCRaw: va.info.EigenPCRaw,
        EigenPCPhred: va.info.EigenPCPhred
    }}''')

vds_1000g_freqs_fields = ['va.freqs.tgp']
vds_hrc_freqs_fields = ['va.freqs.hrc']
vds_gnomad_freqs_fields = ['va.freqs.gnomad']
vds_mgrb_split_freqs_fields = ['va.freqs.mgrb']
vds_dbsnp_fields = ['va.rsid']
vds_clinvar_fields = ['va.clinvar']
vds_cato_fields = ['va.predictions.CATO']
vds_eigen_fields = ['va.predictions.Eigen']
vds_vep_fields = ['va.vep']

# Add successive layers of annotation to the MGRB variant-only vds using the split datasets
(annotate_non_split_from_split(hc, 
                              non_split_vds_path='tmp/05_mgrb_varsonly.vds', 
                              split_vds=vds_1000g_freqs, 
                              annotations=vds_1000g_freqs_fields)
    .write('tmp/05_mgrb_varsonly.1000g.vds')
)

(annotate_non_split_from_split(hc, 
                              non_split_vds_path='tmp/05_mgrb_varsonly.1000g.vds', 
                              split_vds=vds_hrc_freqs, 
                              annotations=vds_hrc_freqs_fields)
    .write('tmp/05_mgrb_varsonly.1000g.hrc.vds')
)

(annotate_non_split_from_split(hc, 
                              non_split_vds_path='tmp/05_mgrb_varsonly.1000g.hrc.vds', 
                              split_vds=vds_gnomad_freqs, 
                              annotations=vds_gnomad_freqs_fields)
    .write('tmp/05_mgrb_varsonly.1000g.hrc.gnomad.vds')
)

(annotate_non_split_from_split(hc, 
                              non_split_vds_path='tmp/05_mgrb_varsonly.1000g.hrc.gnomad.vds', 
                              split_vds=vds_mgrb_split_freqs, 
                              annotations=vds_mgrb_split_freqs_fields)
    .write('tmp/05_mgrb_varsonly.1000g.hrc.gnomad.mgrb_simple.vds')
)

(annotate_non_split_from_split(hc, 
                              non_split_vds_path='tmp/05_mgrb_varsonly.1000g.hrc.gnomad.mgrb_simple.vds', 
                              split_vds=vds_dbsnp, 
                              annotations=vds_dbsnp_fields)
    .write('tmp/05_mgrb_varsonly.1000g.hrc.gnomad.mgrb_simple.dbsnp.vds')
)

(annotate_non_split_from_split(hc, 
                              non_split_vds_path='tmp/05_mgrb_varsonly.1000g.hrc.gnomad.mgrb_simple.dbsnp.vds', 
                              split_vds=vds_clinvar, 
                              annotations=vds_clinvar_fields)
    .write('tmp/05_mgrb_varsonly.1000g.hrc.gnomad.mgrb_simple.dbsnp.clinvar.vds')
)

(annotate_non_split_from_split(hc, 
                              non_split_vds_path='tmp/05_mgrb_varsonly.1000g.hrc.gnomad.mgrb_simple.dbsnp.clinvar.vds', 
                              split_vds=vds_cato, 
                              annotations=vds_cato_fields)
    .write('tmp/05_mgrb_varsonly.1000g.hrc.gnomad.mgrb_simple.dbsnp.clinvar.cato.vds')
)

(annotate_non_split_from_split(hc, 
                              non_split_vds_path='tmp/05_mgrb_varsonly.1000g.hrc.gnomad.mgrb_simple.dbsnp.clinvar.cato.vds', 
                              split_vds=vds_eigen, 
                              annotations=vds_eigen_fields)
    .write('tmp/05_mgrb_varsonly.1000g.hrc.gnomad.mgrb_simple.dbsnp.clinvar.cato.eigen.vds')
)

(annotate_non_split_from_split(hc, 
                              non_split_vds_path='tmp/05_mgrb_varsonly.1000g.hrc.gnomad.mgrb_simple.dbsnp.clinvar.cato.eigen.vds', 
                              split_vds=vds_vep, 
                              annotations=vds_vep_fields)
    .write('tmp/05_mgrb_varsonly.1000g.hrc.gnomad.mgrb_simple.dbsnp.clinvar.cato.eigen.vep.vds')
)


#from hail.representation import Interval
#intervals = [Interval.parse('21:9411850-9411900'), Interval.parse('Y:2649800-2649900'), Interval.parse('X:60000-60050')]

# Perform the final variant-level annotation
#vds_mgrb = vds_mgrb.filter_intervals(intervals)

annotation_vds = hc.read('tmp/05_mgrb_varsonly.1000g.hrc.gnomad.mgrb_simple.dbsnp.clinvar.cato.eigen.vep.vds')
#annotation_vds = annotation_vds.filter_intervals(intervals)
vds_mgrb = vds_mgrb.annotate_variants_vds(annotation_vds, expr='va = vds')

# Recode the variant annotation schema.  Rearrange into the following general structure:
# va: Struct {
#   variant: Struct {
#     Variant-level annotations as structs / vars here
#   },
#   alleles: Array[Struct] {
#     Allele-level annotations as structs / vars here
#   },
#   genotypes: Array[Struct] {
#     Genotype-level annotations as structs / vars here
#   }
# }

# Note the use of va.aIndex and not va.aIndex-1 during the subsetting
# of the va.info.alleleCounts below.  This is deliberate: va.info.alleleCounts
# is derived from oneHotAlleles(v), which includes the reference allele at
# index 0.
vds_mgrb = (vds_mgrb
    .annotate_variants_expr('''
        va.variant = {
            qual: va.qual,
            filters: va.filters,
            metrics: {
                sample_counts: {
                    female:
                        if (v.inYNonPar())
                            0
                        else
                            gs.filter(s => sa.pheno.isFemale).filter(g => g.isCalled()).count(),
                    male:
                        gs.filter(s => !sa.pheno.isFemale).filter(g => g.isCalled()).count() },
                chromatid_counts: {
                    female: va.info.alleleCounts.female.sum(),
                    male: va.info.alleleCounts.male.sum()
                },
                gatk: drop(va.gatk, AC, AF, MLEAC, MLEAF) } },
        va.genotypes = {
            counts: va.info.gtCounts }''')
    .annotate_variants_expr('''
        va.variant.metrics.sample_counts.total = va.variant.metrics.sample_counts.female + va.variant.metrics.sample_counts.male,
        va.variant.metrics.chromatid_counts.total = va.variant.metrics.chromatid_counts.female + va.variant.metrics.chromatid_counts.male''')
    .annotate_alleles_expr('''
        va.alleles = {
            metrics: {
                allele_counts: {
                    female: va.info.alleleCounts.female[va.aIndex],
                    male: va.info.alleleCounts.male[va.aIndex],
                    total: va.info.alleleCounts.total[va.aIndex] },
                allele_frequencies: {
                    female: va.info.alleleCounts.female[va.aIndex] / va.variant.metrics.chromatid_counts.female,
                    male: va.info.alleleCounts.male[va.aIndex] / va.variant.metrics.chromatid_counts.male,
                    total: va.info.alleleCounts.total[va.aIndex] / va.variant.metrics.chromatid_counts.total
                },
                gatk: {
                    AC: va.gatk.AC[va.aIndex-1],
                    AF: va.gatk.AF[va.aIndex-1],
                    MLEAC: va.gatk.AC[va.aIndex-1],
                    MLEAF: va.gatk.AF[va.aIndex-1] } },
            annotation: {
                rsid: va.rsid[va.aIndex-1],
                vep: va.vep[va.aIndex-1],
                clinvar: va.clinvar[va.aIndex-1],
                predictions: {
                    CATO: va.predictions.CATO[va.aIndex-1],
                    Eigen: va.predictions.Eigen[va.aIndex-1] },
                freqs: {
                    tgp: va.freqs.tgp[va.aIndex-1],
                    hrc: va.freqs.hrc[va.aIndex-1],
                    gnomad: va.freqs.gnomad[va.aIndex-1] } } }''')
    .annotate_variants_expr('va = select(va, locus, variant, alleles, genotypes)')
)

# TODO: Add gene-level annotations (eg GAVIN, RVIS)

vds_mgrb.write('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.vds')

