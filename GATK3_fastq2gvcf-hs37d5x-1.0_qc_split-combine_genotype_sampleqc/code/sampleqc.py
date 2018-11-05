#!../../software/pyhail.sh
import pyspark
import hail
from hail.expr import TString, TBoolean

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

import pandas as pd
import pandas.tools.plotting


hc = hail.HailContext(log = 'log/load.log', tmp_dir = 'tmp/hail')

qc_vds = hc.import_vcf('../../GATK3_fastq2gvcf-hs37d5x-1.0_qc_split-combine_genotype/MGRB.phase2qc.hc.norm.vcf.bgz')

sample_table = (hc
    .import_table('../../metadata/MGRB_phase2_metadata.csv', delimiter=',', types={'sampleID': TString(), 'cohort': TString(), 'isFemale': TBoolean()})
    .key_by('sampleID')
)

qc_vds = qc_vds.annotate_samples_table(sample_table, root='sa.pheno')

qc_autosomes_vds = (qc_vds
    .filter_variants_expr('v.isAutosomal()')
    .annotate_samples_expr('sa.hetStats = gs.filter(g => g.isHet() && g.ad[gtk(g.gt)] > 0).map(g => g.ad[gtj(g.gt)] / g.ad[gtk(g.gt)]).stats()')
    .sample_qc()
)
qc_sexcheck_vds = qc_vds.filter_variants_expr('!v.isAutosomal()').filter_multi().variant_qc().filter_variants_expr('va.qc.AF >= 0.05 && va.qc.AF <= 0.95').impute_sex()


qc_autosomes_table = qc_autosomes_vds.samples_table().to_pandas()
qc_sexcheck_table = qc_sexcheck_vds.samples_table().to_pandas()
qc_merged_table = pandas.merge(qc_autosomes_table, qc_sexcheck_table, on='s')
qc_merged_table_forplot = qc_merged_table.loc[:,['sa.qc.callRate', 'sa.qc.dpMean', 'sa.qc.dpStDev', 'sa.qc.gqStDev', 'sa.hetStats.stdev', 'sa.qc.rHetHomVar', 'sa.qc.rTiTv', 'sa.imputesex.Fstat']]
qc_plot = pandas.tools.plotting.scatter_matrix(qc_merged_table_forplot, alpha=0.2, figsize=(16, 16), diagonal='hist')
plt.savefig(r'../qc.png')

qc_merged_table_forexport = qc_merged_table.loc[:,['s', 'sa.qc.callRate', 'sa.qc.dpMean', 'sa.qc.dpStDev', 'sa.qc.gqStDev', 'sa.hetStats.stdev', 'sa.qc.rHetHomVar', 'sa.qc.rTiTv', 'sa.imputesex.Fstat']]
qc_merged_table_forexport.to_csv('../metrics.csv')


# Place samples into tiers as follows:
# 1) Tier 1 -- absolutely top quality sequence.  Very high call rate, very consistent coverage, 
#    excellent genotyping quality, no unusual metrics.
#      callRate > 0.996 &&
#      dpStDev < 8 && 
#      hetStDev < 0.7 && 
#      rHetHomVar < 2 && 
#      (Fstat < 0.2 || Fstat > 0.8) && 
#      (dpStDev < dpMean*0.12 + 2.7) && 
#      (gqStDev < dpMean*-0.68 + 36) && 
#      rTiTv > 3.8 && 
#      rTiTv < 4.3 && 
#      rSingleton < 0.001
# 2) Tier 2 -- probably good sequence.  Fairly high call rate, somewhat consistent coverage,
#    major metrics OK.
#      callRate > 0.98 && 
#      dpStDev < 10 && 
#      hetStDev < 1.0 && 
#      rHetHomVar < 2 && 
#      (Fstat < 0.2 || Fstat > 0.8) && 
#      rSingleton < 0.001
# 3) Tier 3 -- rough but could be made to work.
#      callRate > 0.96 && 
#      hetStDev < 1.0 && 
#      rHetHomVar < 2 &&
#      (Fstat < 0.2 || Fstat > 0.8) && 
#      rSingleton < 0.001
# 4) Tier 4 -- rubbish.  Defined by exclusion: any samples not in one of the above tiers.

good_fstat_samples = set(qc_sexcheck_vds.query_samples('samples.filter(s => sa.imputesex.Fstat < 0.2 || sa.imputesex.Fstat > 0.8).collect()'))
sex_match_samples = set(qc_sexcheck_vds.query_samples('samples.filter(s => sa.imputesex.isFemale == sa.pheno.isFemale).collect()'))
tier_1_samples = set(qc_autosomes_vds.query_samples('''
    samples.filter(s => 
        sa.qc.callRate > 0.996 &&
        sa.qc.dpStDev < 8 &&
        sa.hetStats.stdev < 0.7 &&
        sa.qc.rHetHomVar < 2 &&
        (sa.qc.dpStDev < sa.qc.dpMean*0.12 + 2.7) &&
        (sa.qc.gqStDev < sa.qc.dpMean*(-0.68) + 36) &&
        sa.qc.rTiTv > 3.8 &&
        sa.qc.rTiTv < 4.3 &&
        sa.qc.nSingleton / (sa.qc.nSNP + sa.qc.nInsertion + sa.qc.nDeletion) < 0.001).collect()''')) & good_fstat_samples
tier_2_samples = (set(qc_autosomes_vds.query_samples('''
    samples.filter(s => 
        sa.qc.callRate > 0.98 &&
        sa.qc.dpStDev < 10 &&
        sa.hetStats.stdev < 1.0 &&
        sa.qc.rHetHomVar < 2 &&
        sa.qc.nSingleton / (sa.qc.nSNP + sa.qc.nInsertion + sa.qc.nDeletion) < 0.001).collect()''')) & good_fstat_samples) - tier_1_samples
tier_3_samples = (set(qc_autosomes_vds.query_samples('''
    samples.filter(s => 
        sa.qc.callRate > 0.96 &&
        sa.hetStats.stdev < 1.0 &&
        sa.qc.rHetHomVar < 2 &&
        sa.qc.nSingleton / (sa.qc.nSNP + sa.qc.nInsertion + sa.qc.nDeletion) < 0.001).collect()''')) & good_fstat_samples) - tier_1_samples - tier_2_samples
tier_4_samples = set(qc_autosomes_vds.query_samples('samples.collect()')) - tier_1_samples - tier_2_samples - tier_3_samples


tier_1_sexmatch_samples = tier_1_samples & sex_match_samples
tier_2_sexmatch_samples = tier_2_samples & sex_match_samples
tier_3_sexmatch_samples = tier_3_samples & sex_match_samples
tier_1_sexmismatch_samples = tier_1_samples - tier_1_sexmatch_samples
tier_2_sexmismatch_samples = tier_2_samples - tier_2_sexmatch_samples
tier_3_sexmismatch_samples = tier_3_samples - tier_3_sexmatch_samples

open('../tier1.sexmatch.sample_list', 'wt').write('\n'.join(sorted(tier_1_sexmatch_samples)) + '\n')
open('../tier2.sexmatch.sample_list', 'wt').write('\n'.join(sorted(tier_2_sexmatch_samples)) + '\n')
open('../tier3.sexmatch.sample_list', 'wt').write('\n'.join(sorted(tier_3_sexmatch_samples)) + '\n')
open('../tier1.sexmismatch.sample_list', 'wt').write('\n'.join(sorted(tier_1_sexmismatch_samples)) + '\n')
open('../tier2.sexmismatch.sample_list', 'wt').write('\n'.join(sorted(tier_2_sexmismatch_samples)) + '\n')
open('../tier3.sexmismatch.sample_list', 'wt').write('\n'.join(sorted(tier_3_sexmismatch_samples)) + '\n')
open('../tier4.all.sample_list', 'wt').write('\n'.join(sorted(tier_4_samples)) + '\n')
