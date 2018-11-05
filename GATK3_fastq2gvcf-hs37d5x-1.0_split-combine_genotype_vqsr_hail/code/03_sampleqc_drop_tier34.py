#!./bin/pyhail.sh
import pyspark
import hail

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

import pandas as pd
import pandas.tools.plotting


hc = hail.HailContext(log = 'log/03_sampleqc_drop_tier34.log', tmp_dir = 'tmp/hail')

vds = hc.read('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.vds')

qcregions = hail.KeyTable.import_bed('../../locus-annotations/regions/qc.bed.gz')

# Simple sample QC, following the procedure used for earlier targeted QC (based on QC SNP
# chip loci), but expanded to whole-genome callable regions from locus-annotations_hail/regions/qc.bed.gz.
# Do not perform sex imputation and checking again, as this is expected to be good enough 
# based on the initial 'SNP-chip' QC.
# 
# Some thresholds have been changed from the earlier targeted QC, to reflect the
# broader selection of loci in this pass:
# 1) Tier 1 -- absolutely top quality sequence.  Very high call rate, very consistent coverage, 
#    excellent genotyping quality, no unusual metrics.
#      This                                      SNP QC
#      callRate > 0.996 &&                       Same
#      dpStDev < 6 &&                            dpStDev < 8
#      hetStDev < 0.7 &&                         Same
#      rHetHomVar < 2 &&                         Same
#      (Not tested)                              (Fstat < 0.2 || Fstat > 0.8)
#      (dpStDev < dpMean*-0.05 + 7.5)            (dpStDev < dpMean*0.12 + 2.7)
#      (gqStDev < dpMean*-0.68 + 38) &&          (gqStDev < dpMean*-0.68 + 36)
#      rTiTv > 2.05 &&                           rTiTv > 3.8
#      rTiTv < 2.15 &&                           rTiTv < 4.3
#      (Not tested)                              rSingleton < 0.001
# 2) Tier 2 -- probably good sequence.  Fairly high call rate, somewhat consistent coverage,
#    major metrics OK.
#      This                                      SNP QC
#      callRate > 0.98 &&                        Same
#      dpStDev < 10 &&                           Same
#      hetStDev < 1.0 &&                         Same
#      rHetHomVar < 2 &&                         Same
#      (Not tested)                              (Fstat < 0.2 || Fstat > 0.8)
#      (Not tested)                              rSingleton < 0.001
# 3) Tier 3 -- rough but could be made to work.
#      callRate > 0.96 &&                        Same
#      hetStDev < 1.0 &&                         Same
#      rHetHomVar < 2 &&                         Same
#      (Not tested)                              (Fstat < 0.2 || Fstat > 0.8)
#      (Not tested)                              rSingleton < 0.001
# 4) Tier 4 -- rubbish.  Defined by exclusion: any samples not in one of the above tiers.
vds_qc = (vds
    .filter_variants_table(qcregions)
    .sample_qc()
    .annotate_samples_expr('sa.qc.hetStDev = gs.filter(g => g.isHet() && g.ad[gtk(g.gt)] > 0).map(g => g.ad[gtj(g.gt)] / g.ad[gtk(g.gt)]).stats().stdev')
)

vds_qc = vds_qc.annotate_samples_expr('sa.qc.hetStDev = gs.filter(g => g.isHet() && g.ad[gtk(g.gt)] > 0).map(g => g.ad[gtj(g.gt)] / g.ad[gtk(g.gt)]).stats().stdev')

qc_table = vds_qc.samples_table().to_pandas()
qc_table_forplot = qc_table.loc[:,['sa.qc.callRate', 'sa.qc.dpMean', 'sa.qc.dpStDev', 'sa.qc.gqStDev', 'sa.qc.hetStDev', 'sa.qc.rHetHomVar', 'sa.qc.rTiTv']]
qc_plot = pandas.tools.plotting.scatter_matrix(qc_table_forplot, alpha=0.2, figsize=(16, 16), diagonal='hist')
plt.savefig(r'../qc.svg', format="svg")

qc_table_forexport = qc_table.loc[:,['s', 'sa.qc.callRate', 'sa.qc.dpMean', 'sa.qc.dpStDev', 'sa.qc.gqStDev', 'sa.qc.hetStDev', 'sa.qc.rHetHomVar', 'sa.qc.rTiTv']]
qc_table_forexport.to_csv('../03_metrics.csv')

tier_1_samples = set(vds_qc.query_samples('''
    samples.filter(s => 
        sa.qc.callRate > 0.996 &&
        sa.qc.dpStDev < 6 &&
        sa.qc.hetStDev < 0.7 &&
        sa.qc.rHetHomVar < 2 &&
        (sa.qc.dpStDev < sa.qc.dpMean*(-0.05) + 7.5) &&
        (sa.qc.gqStDev < sa.qc.dpMean*(-0.68) + 38) &&
        sa.qc.rTiTv > 2.05 &&
        sa.qc.rTiTv < 2.15).collect()'''))
tier_2_samples = (set(vds_qc.query_samples('''
    samples.filter(s => 
        sa.qc.callRate > 0.98 &&
        sa.qc.dpStDev < 10 &&
        sa.qc.hetStDev < 1.0 &&
        sa.qc.rHetHomVar < 2).collect()'''))) - tier_1_samples
tier_3_samples = (set(vds_qc.query_samples('''
    samples.filter(s => 
        sa.qc.callRate > 0.96 &&
        sa.qc.hetStDev < 1.0 &&
        sa.qc.rHetHomVar < 2).collect()'''))) - tier_1_samples - tier_2_samples
tier_4_samples = set(vds_qc.query_samples('samples.collect()')) - tier_1_samples - tier_2_samples - tier_3_samples



vds = hc.read('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.vds')

open('../tier1.sample_list', 'wt').write('\n'.join(sorted(tier_1_samples)) + '\n')
open('../tier2.sample_list', 'wt').write('\n'.join(sorted(tier_2_samples)) + '\n')
open('../tier3.sample_list', 'wt').write('\n'.join(sorted(tier_3_samples)) + '\n')
open('../tier4.sample_list', 'wt').write('\n'.join(sorted(tier_4_samples)) + '\n')

# Filter the original VDS: Exclude tier 3 and 4 samples, and mark remaining samples as either tier 1 or 2.
# Also filter alleles and variants to remove entries that no longer are supported by any samples.
tier1_table = hc.import_table('../tier1.sample_list', no_header=True).annotate('sampleID = f0').key_by('sampleID')
#tier2_table = hc.import_table('../tier2.sample_list', no_header=True).annotate('sampleID = f0').key_by('sampleID')

#    .annotate_variants_expr('''va.info = {
#        callRate: gs.fraction(g => g.isCalled),
#        isSingleton: gs.filter(g => g.isCalledNonRef).count() == 1,
#        NS: gs.filter(g => g.isCalled).count()}''')
#    .annotate_alleles_expr('va.info.AC = gs.filter(g => g.isHet).count() + gs.filter(g => g.isHomVar).count()*2')
#    .annotate_variants_expr('va.info.AF = va.info.AC / (va.info.NS*2)')
#    .filter_alleles('va.info.AC[aIndex - 1] == 0', annotation='''
#        va.info.AC = aIndices[1:].map(i => va.info.AC[i - 1]),
#        va.info.AF = aIndices[1:].map(i => va.info.AF[i - 1])''',
#        keep=False, subset=True, keep_star=False)

vds = (vds
    .filter_samples_list(list(tier_1_samples) + list(tier_2_samples))
    .filter_variants_expr('gs.filter(g => g.isCalledNonRef).count() > 0')
    .naive_coalesce(5000)
    .annotate_samples_table(tier1_table, root='sa.tier1')
    .annotate_samples_expr('sa.tier = if (isMissing(sa.tier1)) 2 else 1')
    .annotate_samples_expr('sa = drop(sa, tier1)')
    .write('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.vds')
)
