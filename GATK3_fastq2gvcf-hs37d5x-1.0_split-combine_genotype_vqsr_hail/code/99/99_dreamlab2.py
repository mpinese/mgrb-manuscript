#!../../software/pyhail.sh
import hail
from hail.representation import Interval
from hail.expr import TString, TBoolean, TFloat, TInt


hc = hail.HailContext(log = 'log/99_dreamlab2.log', tmp_dir = 'tmp/hail')

vds = hc.read('../MGRB.phase2.tier12.match.vqsr.minrep.vds')

# Chr22 only, rough variant quality filters
vds = (vds
    .filter_intervals(Interval.parse('22'), keep=True)
    .filter_variants_expr('va.filters.isEmpty()', keep=True)
    .split_multi()
    .variant_qc()
    .filter_variants_expr('''
        v.altAllele.isSNP &&
        va.qc.callRate >= 0.99 && 
        va.qc.dpMean >= 20 && va.qc.dpMean <= 60 && 
        va.qc.dpStDev < 8 && 
        va.filters.isEmpty() && 
        va.qc.AF >= 0.05 && va.qc.AF <= 0.95''')
)

# Drop samples with poor metrics on these filtered variants.
vds = (vds
    .sample_qc()
    .filter_samples_expr('sa.qc.callRate >= 0.985')
)

# Drop samples with no phenotype for SBPMean, HtMtrs, WtKgs, AMD, or GlcmmolL
vds = (vds
    .filter_samples_expr('''
        isnan(sa.pheno.SBPMean) ||
        isMissing(sa.pheno.SBPMean) ||
        isnan(sa.pheno.HtMtrs) ||
        isMissing(sa.pheno.HtMtrs) ||
        isnan(sa.pheno.WtKgs) ||
        isMissing(sa.pheno.WtKgs) ||
        isMissing(sa.pheno.AMD) ||
        isnan(sa.pheno.GlcmmolL) ||
        isMissing(sa.pheno.GlcmmolL)''', keep=False)
    .annotate_samples_expr('sa.pheno = select(sa.pheno, SBPMean, HtMtrs, WtKgs, AMD, GlcmmolL)')
)


vds = vds.repartition(280)

# PCs for covariates
vds = vds.pca('sa.scores', k=5)

# Export the allele count table
vds.make_table('v = v', ['aac = if (g.isHet()) 1 else if (g.isHomVar()) 2 else 0']).export('../MGRB.phase2.tier12.match.vqsr.minrep.dreamlab2chr22.geno.tsv')

# Export the covariates
vds.export_samples('../MGRB.phase2.tier12.match.vqsr.minrep.dreamlab2chr22.covar.tsv', 'sample = s, sa.scores.*')

# Export the phenotypes
vds.export_samples('../MGRB.phase2.tier12.match.vqsr.minrep.dreamlab2chr22.pheno.tsv', 'sample = s, sa.pheno.*')

