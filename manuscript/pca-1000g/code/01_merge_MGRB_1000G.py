#!./bin/pyhail.sh
import hail
import os
import pprint
from hail.expr import TVariant


hc = hail.HailContext(log = 'log/01_merge_mgrb_1000G.log', tmp_dir = 'tmp/hail')

vds_mgrb = hc.read('../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.vds')

vds_1000g_highqual = hc.read('../../../databases/1000G_hail/1000G.minrep.qcregions.common.dbsnp.biallelic.hwe.vds')



# Merge MGRB and 1000 genomes at the genotype level to create a dataset for
# population comparisons.  Due to the filtering used to create vds_1000g_highqual,
# this will be restricted to unambiguous autosomal SNPs with RSIDs only,
# with no evidence of HWE deviation within 1000G homogeneous cohorts.

# Filter MGRB to high quality SNVs, and split multi-allelics
vds_mgrb_highqual = (vds_mgrb
    .hardcalls()
    .split_multi()
    .filter_variants_expr('v.altAllele().isSNP() && va.variant.filters.isEmpty() && va.locus.tier != 3')
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

# Perform variant QC and drop variants with very high or low frequencies.
vds_mgrb_1000g = vds_mgrb_1000g.variant_qc().filter_variants_expr('va.qc.AF >= 0.05 && va.qc.AF <= 0.95')

# Save as an intermediate merged dataset, without pruning
vds_mgrb_1000g.write('../01_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.1000G_commonSNVs.vds')

# Also create an LD-pruned version of the MGRB-1000G cohort.
# Use a prespecified LD pruned variant list, if present.  If not, create one.
if os.path.isfile('../01_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.1000G_commonSNVs.ldpruned.variant_list'):
    print('Pruning with prespecified LD loci...')
    variants_kt = (hc
        .import_table('../01_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.1000G_commonSNVs.ldpruned.variant_list', no_header=True, types={'f0': TVariant()})
        .key_by('f0'))
    vds_mgrb_1000g_pruned = vds_mgrb_1000g.filter_variants_table(variants_kt)
else:
    print('Finding and pruning LD loci...')
    vds_mgrb_1000g_pruned = vds_mgrb_1000g.ld_prune(r2 = 0.1)
    vds_mgrb_1000g_pruned.export_variants('../01_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.1000G_commonSNVs.ldpruned.variant_list', 'v')

vds_mgrb_1000g_pruned = vds_mgrb_1000g_pruned.repartition(200)

vds_mgrb_1000g_pruned.write('../01_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.1000G_commonSNVs.ldpruned.vds')


