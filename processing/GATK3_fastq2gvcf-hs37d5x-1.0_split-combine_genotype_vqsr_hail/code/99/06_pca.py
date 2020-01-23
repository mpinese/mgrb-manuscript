#!./bin/pyhail.sh
import hail
import os
import pprint
from hail.expr import TVariant



hc = hail.HailContext(log = 'log/05_merge_with_ref_cohorts.log', tmp_dir = 'tmp/hail')

vds_mgrb = hc.read('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.vds')

"""
vds_1000g_highqual = hc.read('../../1000G_hail/1000G.minrep.qcregions.common.dbsnp.biallelic.hwe.vds')



# Merge MGRB and 1000 genomes at the genotype level to create a dataset for
# population comparisons.  Due to the filtering used to create vds_1000g_highqual,
# this will be restricted to unambiguous autosomal SNPs with RSIDs only, that
# with no evidence of HWE deviation within 1000G homogeneous cohorts.

# Filter MGRB to high quality SNVs, and split multi-allelics
vds_mgrb_highqual = (vds_mgrb
    .hardcalls()
    .split_multi()
    .filter_variants_expr('v.altAllele().isSNP() && va.filters.isEmpty() && va.locus.tier != 3')
    .naive_coalesce(1000)
)

# Prior to the join we must bring sample annotations into alignment
vds_1000g_highqual = vds_1000g_highqual.annotate_samples_expr('''
    sa.pheno = {
        manifest: NA: String,
        FRID: NA: String,
        externalID: NA: String,
        cohort: sa.pheno.cohort,
        YOB: NA: Int,
        SBPMean: NA: Int,
        HtMtrs: NA: Float,
        WtKgs: NA: Float,
        AbdoCircCms: NA: Int,
        GlcmmolL: NA: Float,
        AMD: NA: Boolean,
        isFemale: sa.pheno.isFemale,
        treatedForHighBP: NA: Boolean,
        treatedForHighChol: NA: Boolean,
        population: sa.pheno.population,
        superPopulation: sa.pheno.superPopulation
    },
    sa.qc = {tier: NA: Int}''')

vds_mgrb_highqual = (vds_mgrb_highqual.annotate_samples_expr('''
    sa.pheno.population = "MGRB",
    sa.pheno.superPopulation = "MGRB"''')
)

# Merge MGRB with 1000G
vds_mgrb_1000g = vds_mgrb_highqual.join(vds_1000g_highqual).repartition(1000)

# Save as an intermediate merged dataset, without pruning
vds_mgrb_1000g.write('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.1000G_commonSNVs.vds')

# Also create an LD-pruned version of the MGRB-1000G cohort.
# Use a prespecified LD pruned variant list, if present.  If not, create one.
if os.path.isfile('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.1000G_commonSNVs.ldpruned.variant_list'):
    print('Using existing LD pruning loci...')
    variants_kt = (hc
        .import_table('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.1000G_commonSNVs.ldpruned.variant_list', no_header=True, types={'f0': TVariant()})
        .key_by('f0'))
    vds_mgrb_1000g_pruned = vds_mgrb_1000g.filter_variants_table(variants_kt)
else:
    print('Finding LD pruning loci...')
    vds_mgrb_1000g_pruned = vds_mgrb_1000g.ld_prune(r2 = 0.1)
    vds_mgrb_1000g_pruned.export_variants('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.1000G_commonSNVs.ldpruned.variant_list', 'v')

vds_mgrb_1000g_pruned = vds_mgrb_1000g_pruned.repartition(200)

vds_mgrb_1000g_pruned.write('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.1000G_commonSNVs.ldpruned.vds')






# Generate a temporary frequency-only version of the MGRB.  This is in part to
# act as a QC for the autosomal MGRB frequencies already calculated, as well as
# to provide sex naive frequency estimates for better comparison to other cohorts
# which don't do sex-aware AF calculation.
vss_mgrb.split_multi().annotate_variants_expr('va = {}').variant_qc().drop_samples().write('tmp/05_mgrb_simpleafs.vds')





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

    variant_annotated_vds = (
        hc.read(non_split_vds_path, drop_samples=True)
        .annotate_variants_expr('va.variant = v')
        .split_multi()
    )

    ann_agg_codes = ["`%s` = index(va.map(x => {val: %s, aIndex: va.aIndex}).collect(), aIndex)" % (a, a) for a in annotations]
    agg = (
        split_vds
            .annotate_variants_vds(variant_annotated_vds, 'va.variant = vds.variant, va.aIndex = vds.aIndex')
            .filter_variants_expr('isDefined(va.variant)')
            .variants_table()
            .aggregate_by_key('variant = va.variant', ann_agg_codes + ['nAltAlleles = va.map(x => x.aIndex).max()'])
    )

    ann_codes = ['%s = let x = table.`%s` in' \
                 ' range(table.nAltAlleles).map(i => if(x.contains(i+1)) x[i+1].val else NA: %s)' % (a, a, b)
                 for (a, b) in zip(annotations, ann_types)]

    return (
        hc.read(non_split_vds_path)
            .annotate_variants_table(agg, expr=",".join(ann_codes))
    )


# Create the MGRB annotation-only VDS starting point
#######################################################################################################################################################################
#vds_mgrb.drop_samples().naive_coalesce(200).write('tmp/05_mgrb_varsonly.vds')

# Load the source allele frequency databases
#vds_1000g_freqs = hc.read('../../1000G_hail/1000G.split.minrep.freqs.vds')
#vds_hrc_freqs = hc.read('../../HRC_hail/HRC.minrep.vds')
#vds_gnomad_freqs = hc.read('../../GnomAD_hail/gnomad.genomes.r2.0.1.sites.combined.split.minrep.vds')
vds_mgrb_split_freqs = hc.read('tmp/05_mgrb_simpleafs.vds')

# Relabel the annotations in the source databases to avoid collisions

vds_1000g_freqs = vds_1000g_freqs.annotate_variants_expr('''
    va.freqs = {tgp: {
        AC: va.info.AC[va.aIndex - 1],
        AF: va.info.AF[va.aIndex - 1],
        NS: va.info.NS,
        AN: va.info.AN,
        AF_EAS: va.info.EAS_AF[va.aIndex - 1],
        AF_EUR: va.info.EUR_AF[va.aIndex - 1],
        AF_AFR: va.info.AFR_AF[va.aIndex - 1],
        AF_AMR: va.info.AMR_AF[va.aIndex - 1],
        AF_SAS: va.info.SAS_AF[va.aIndex - 1],
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
vds_gnomad_freqs = vds_gnomad_freqs.annotate_variants_expr('va.freqs = {gnomad: va.gnomad}')

vds_mgrb_split_freqs = vds_mgrb_split_freqs.annotate_variants_expr('''
    va.info = {
        AC_Hail: va.qc.AC,
        AF_Hail: va.qc.AF,
        NS_Hail: va.qc.nCalled
    }''')

vds_1000g_freqs_fields = ['va.freqs.tgp.' + s for s in ['AC', 'AF', 'NS', 'AN', 'AF_EAS', 'AF_EUR', 'AF_AFR', 'AF_AMR', 'AF_SAS', 'AC_Hail', 'AF_Hail', 'NS_Hail']]
vds_hrc_freqs_fields = ['va.freqs.hrc.' + s for s in ['AC', 'AN', 'AF', 'AC_EXCLUDING_1000G', 'AN_EXCLUDING_1000G', 'AF_EXCLUDING_1000G']]
vds_gnomad_freqs_fields = ['va.freqs.gnomad.' + s for s in ['AC', 'AF', 'AN', 'AC_AFR', 'AC_AMR', 'AC_ASJ', 'AC_EAS', 'AC_FIN', 'AC_NFE', 'AC_OTH', 'AN_AFR', 'AN_AMR', 'AN_ASJ', 'AN_EAS', 'AN_FIN', 'AN_NFE', 'AN_OTH', 'AF_AFR', 'AF_AMR', 'AF_ASJ', 'AF_EAS', 'AF_FIN', 'AF_NFE', 'AF_OTH']]
vds_mgrb_split_freqs_fields = ['va.info.' + s for s in ['AC_Hail', 'AF_Hail', 'NS_Hail']]

# Add successive layers of annotation to this vds using the split datasets

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
"""

# TODO: Add ClinVar, dbSNP rsIDs.

# Also TODO: Split off all this annotation code to a separate file.  Place
# the 1000 G merge afterwards.  Consider VEP in-between.  Would be nice to
# have the annotations locked down before subsetting.

# Perform the final annotation
(vds_mgrb
    .annotate_variants_vds('tmp/05_mgrb_varsonly.1000g.hrc.gnomad.mgrb_simple.vds', expr='va = vds.va')
    .write('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.popfreqs.vds')
)
