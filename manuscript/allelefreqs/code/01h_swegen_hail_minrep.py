#!./bin/pyhail.sh
import hail

hc = hail.HailContext(log = 'log/01h_swegen_hail_minrep.log', tmp_dir = 'tmp/hail')

vds = hc.import_vcf('tmp/swegen_split.vcf.bgz')

(vds
    .filter_variants_expr('"^[0-9]+$" ~ v.contig')
    .min_rep()
    .export_variants('tmp/swegen_split_minrep_autosomes.tsv', 'v.contig, v.start, v.ref, v.alt(), va.info.nRR, va.info.nRA, va.info.nAA, va.info.nmissing')
)
