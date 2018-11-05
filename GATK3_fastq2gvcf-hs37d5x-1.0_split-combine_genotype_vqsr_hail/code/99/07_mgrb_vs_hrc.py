#!./bin/pyhail.sh
import hail
import subprocess

hc = hail.HailContext(log = 'log/07_mgrb_vs_hrc.log', tmp_dir = 'tmp/hail')

vds_mgrb = hc.read('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.vds')
vds_hrc = hc.read('../../HRC_hail/HRC.minrep.vds')

# Filter MGRB to high quality SNVs, and split multi-allelics
vds_mgrb_hrc = (vds_mgrb
    .hardcalls()
    .split_multi()
    .filter_variants_expr('v.altAllele().isSNP()')
    .variant_qc()
    .annotate_variants_vds(vds_hrc, expr = '''
        va.rsid = vds.rsid,
        va.hrc = vds.info''')
)

vds_mgrb_hrc.export_variants('../mgrb_hrc_variant_freqs.tsv', 'variant = v, rsid = va.rsid, filters = va.filters.mkString(","), tier = va.locus.tier, mgrb_AF = va.qc.AF, HRC_AF = va.hrc.AF, HRC_AF_EXCLUDING_1000G = va.hrc.AF_EXCLUDING_1000G')

subprocess.call('R --slave -e \'freqs = read.table("../mgrb_hrc_variant_freqs.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE); saveRDS(freqs, "../mgrb_hrc_variant_freqs.rds")\'', shell=True)
