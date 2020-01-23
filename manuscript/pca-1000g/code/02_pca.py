#!./bin/pyhail.sh
import hail
import pprint

hc = hail.HailContext(log = 'log/02_pca.log', tmp_dir = 'tmp/hail')

vds_mgrb_1000g_pruned = hc.read('../01_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.1000G_commonSNVs.ldpruned.vds')

vds_mgrb_1000g_pruned_pca = vds_mgrb_1000g_pruned.pca('sa.scores', 'va.loadings', 'global.evals', 20)

vds_mgrb_1000g_pruned_pca.write('../02_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.1000G_commonSNVs.ldpruned.pca.vds')

vds_mgrb_1000g_pruned_pca.export_samples('../02_MGRB_1000G_pca.scores.tsv', 'sample = s, sa.pheno.*, sa.qc.*, sa.scores.*')
vds_mgrb_1000g_pruned_pca.export_variants('../02_MGRB_1000G_pca.loadings.tsv', 'variant = v, va.loadings.*')
open('../02_MGRB_1000G_pca.evals.txt', 'wt').write(pprint.pformat(vds_mgrb_1000g_pruned_pca.globals.evals))
