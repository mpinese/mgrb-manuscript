#!./bin/pyhail.sh
import hail
from hail.expr import TString, TBoolean, TFloat, TInt

hc = hail.HailContext(log = 'log/03_import_highqual_common_to_hail.log', tmp_dir = 'tmp/hail')

vds = hc.import_vcf('../phase3.normed.vcf/*.vcf.bgz')

qcregions = hail.KeyTable.import_bed('../../../locus-annotations/regions/qc.bed.gz')

sample_table = (hc
    .import_table('../phase3.vcf/integrated_call_samples_v3.20130502.ALL.panel', delimiter='\t', types={
        'sample': TString(),
        'pop': TString(),
        'super_pop':TString(),
        'gender':TString(),
    })
    .annotate('isFemale = gender == "female"')
    .annotate('cohort = "1000G"')
    .rename({'sample': 'sampleID', 'pop': 'population', 'super_pop': 'superPopulation'})
    .select(['sampleID', 'cohort', 'population', 'superPopulation', 'isFemale'])
    .key_by('sampleID')
)

# Split and run variant_qc to get allele frequencies
vds_split = vds.split_multi().min_rep()
print(vds_split.summarize())
vds_split = vds_split.variant_qc().drop_samples()
print(vds_split.summarize())
vds_split = vds_split.annotate_variants_expr('''
    va.info.CIEND = va.info.CIEND[va.aIndex-1],
    va.info.CIPOS = va.info.CIPOS[va.aIndex-1],
    va.info.MC = va.info.MC[va.aIndex-1],
    va.info.MEINFO = va.info.MEINFO[va.aIndex-1],
    va.info.SVLEN = va.info.SVLEN[va.aIndex-1],
    va.info.AC = va.info.AC[va.aIndex-1],
    va.info.AF = va.info.AF[va.aIndex-1],
    va.info.EAS_AF = va.info.EAS_AF[va.aIndex-1],
    va.info.EUR_AF = va.info.EUR_AF[va.aIndex-1],
    va.info.AFR_AF = va.info.AFR_AF[va.aIndex-1],
    va.info.AMR_AF = va.info.AMR_AF[va.aIndex-1],
    va.info.SAS_AF = va.info.SAS_AF[va.aIndex-1]''')
vds_split.repartition(1000).write('../1000G.split.minrep.freqs.vds')

# Generate a very high quality variant set with homogeneous structure 
# within populations, for PCA against MGRB.
# The alt() and ref comparisons below remove strand-ambiguous SNPs.
vds = (vds
    .filter_variants_table(qcregions)
    .filter_multi()
    .min_rep()
    .annotate_samples_table(sample_table, root='sa.pheno')
    .filter_variants_expr('''
        v.isAutosomal &&
        va.filters.isEmpty() &&
        va.info.VT.forall(v => v == "SNP") &&
        va.info.AF.forall(v => v >= 0.05) &&
        va.info.AF.forall(v => v <= 0.95) &&
        va.rsid != "." &&
        !((v.alt() == "G" && v.ref == "C") ||
          (v.alt() == "C" && v.ref == "G") ||
          (v.alt() == "T" && v.ref == "A") ||
          (v.alt() == "A" && v.ref == "T")) &&
        (v.ref == "A" || v.ref == "C" || v.ref == "G" || v.ref == "T") &&
        (v.alt() == "A" || v.alt() == "C" || v.alt() == "G" || v.alt() == "T") &&
        gs.filter(g => sa.pheno.population == "MSL").hardyWeinberg().pHWE > 1e-2/17 &&
        gs.filter(g => sa.pheno.population == "BEB").hardyWeinberg().pHWE > 1e-2/17 &&
        gs.filter(g => sa.pheno.population == "GBR").hardyWeinberg().pHWE > 1e-2/17 &&
        gs.filter(g => sa.pheno.population == "CDX").hardyWeinberg().pHWE > 1e-2/17 &&
        gs.filter(g => sa.pheno.population == "CEU").hardyWeinberg().pHWE > 1e-2/17 &&
        gs.filter(g => sa.pheno.population == "FIN").hardyWeinberg().pHWE > 1e-2/17 &&
        gs.filter(g => sa.pheno.population == "KHV").hardyWeinberg().pHWE > 1e-2/17 &&
        gs.filter(g => sa.pheno.population == "LWK").hardyWeinberg().pHWE > 1e-2/17 &&
        gs.filter(g => sa.pheno.population == "ITU").hardyWeinberg().pHWE > 1e-2/17 &&
        gs.filter(g => sa.pheno.population == "STU").hardyWeinberg().pHWE > 1e-2/17 &&
        gs.filter(g => sa.pheno.population == "CHB").hardyWeinberg().pHWE > 1e-2/17 &&
        gs.filter(g => sa.pheno.population == "JPT").hardyWeinberg().pHWE > 1e-2/17 &&
        gs.filter(g => sa.pheno.population == "CHS").hardyWeinberg().pHWE > 1e-2/17 &&
        gs.filter(g => sa.pheno.population == "IBS").hardyWeinberg().pHWE > 1e-2/17 &&
        gs.filter(g => sa.pheno.population == "TSI").hardyWeinberg().pHWE > 1e-2/17 &&
        gs.filter(g => sa.pheno.population == "YRI").hardyWeinberg().pHWE > 1e-2/17 &&
        gs.filter(g => sa.pheno.population == "GWD").hardyWeinberg().pHWE > 1e-2/17''')
    .annotate_variants_expr('''va.info = {
        AC: va.info.AC[0],
        AF: va.info.AF[0],
        NS: va.info.NS,
        AN: va.info.AN,
        EAS_AF: va.info.EAS_AF[0],
        EUR_AF: va.info.EUR_AF[0],
        AFR_AF: va.info.AFR_AF[0],
        AMR_AF: va.info.AMR_AF[0],
        SAS_AF: va.info.SAS_AF[0]}''')
    .repartition(1000).write('../1000G.minrep.qcregions.common.dbsnp.biallelic.hwe.vds')
)
