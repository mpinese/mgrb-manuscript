#!./bin/pyhail.sh
import pyspark
import hail
from hail.representation import Interval
from hail import KeyTable
import os.path

hc = hail.HailContext(log = 'log/07_plink_export_45andup_subpops.log', tmp_dir = 'tmp/hail')

interval_list = [Interval.parse('%d' % chrom) for chrom in range(1, 22)]
tier1_bed = KeyTable.import_bed('../../locus-annotations/source_data/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed')

(hc
    .read('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.vds')
    .filter_intervals(interval_list, keep=True)
    .annotate_variants_table(tier1_bed, root='va.tier1bed')
    .filter_variants_expr('va.tier1bed == true', keep=True)
    .annotate_variants_expr('va = {}')
    .filter_samples_expr('"^B" ~ s')
    .split_multi()
    .min_rep()
    .write('tmp/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.45andUp.GiaB_HCR.noannot.split.minrep.vds')
)

sample_lists = [
    '../45andup_followup_qcpass_anycancer_mf.sample_list',
    '../45andup_followup_qcpass_breastcancer_f.sample_list',
    '../45andup_followup_qcpass_breastcancer_mf.sample_list',
    '../45andup_followup_qcpass_colorectalcancer_mf.sample_list',
    '../45andup_followup_qcpass_melanomacancer_mf.sample_list',
    '../45andup_followup_qcpass_nocancer_f.sample_list',
    '../45andup_followup_qcpass_nocancer_mf.sample_list',
    '../45andup_followup_qcpass_nocancer_m.sample_list',
    '../45andup_followup_qcpass_nonmelskincancer_mf.sample_list',
    '../45andup_followup_qcpass_prostatecancer_m.sample_list']

for sample_list in sample_lists:
    sample_list_name = os.path.splitext(os.path.basename(sample_list))[0]
    print 'Sample %s' % sample_list_name
    samples = [s.strip() for s in open(sample_list)]
    (hc
        .read('tmp/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.45andUp.GiaB_HCR.noannot.split.minrep.vds')
        .filter_samples_list(samples)
        .export_plink('../plink/07_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.%s.GiaB_HCR.split.minrep' % sample_list_name)
    )
