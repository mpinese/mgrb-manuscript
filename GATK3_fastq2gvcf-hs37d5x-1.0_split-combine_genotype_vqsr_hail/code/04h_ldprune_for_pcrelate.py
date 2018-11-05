#!./bin/pyhail.sh
import hail
import os
import subprocess
from hail.expr import TVariant


hc = hail.HailContext(log = 'log/04h_ldprune_for_pcrelate.log', tmp_dir = 'tmp/hail')

vds_mgrb = hc.read('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.vds')

# Filter MGRB to high quality SNVs, and remove multi-allelics and rare variants
vds_mgrb = (vds_mgrb
    .hardcalls()
    .filter_multi()
    .filter_variants_expr('v.altAllele().isSNP() && va.filters.isEmpty() && va.locus.tier != 3 && v.isAutosomal()')
    .variant_qc()
    .filter_variants_expr('va.qc.AF >= 0.05 && va.qc.AF <= 0.95')
)
print(vds_mgrb.summarize())

# LD prune and save
# Use a prespecified LD pruned variant list, if present.  If not, create one.
if os.path.isfile('../ldpruned/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.ldpruned.variant_list'):
    print('Using existing LD pruning loci...')
    variants_kt = (hc
        .import_table('../ldpruned/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.ldpruned.variant_list', no_header=True, types={'f0': TVariant()})
        .key_by('f0'))
    vds_mgrb_pruned = vds_mgrb.filter_variants_table(variants_kt)
else:
    print('Finding LD pruning loci...')
    vds_mgrb_pruned = vds_mgrb.ld_prune(r2=0.1)
    vds_mgrb_pruned.export_variants('../ldpruned/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.ldpruned.variant_list', 'v')

vds_mgrb_pruned = vds_mgrb_pruned.repartition(200)
print(vds_mgrb_pruned.summarize())

#vds_mgrb_pruned.export_plink('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.ldpruned.plink')
# Export as a VCF, then get PLINK to do the VCF -> BED conversion.
# Reason is that Hail chooses to use its chr:pos:ref:alt convention
# to name alleles, and this causes problems to downstream tools.
# The vds -> vcf -> bed workflow circumvents this issue.
vds_mgrb_pruned.annotate_variants_expr('va = {rsid: va.rsid}').export_vcf('../ldpruned/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.ldpruned.vcf.bgz')
subprocess.call("./bin/plink --vcf ../ldpruned/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.ldpruned.vcf.bgz --out ../ldpruned/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.ldpruned.plink", shell=True)
