#!./bin/pyhail.sh
import hail

hc = hail.HailContext(log = 'log/05h_gnomad2af_tbl.log', tmp_dir = 'tmp/hail')

vds_gnomad_autosomes = hc.read('../gnomad.genomes.r2.0.1.sites.autosomes.vds')

(vds_gnomad_autosomes
    .annotate_variants_expr('va = select(va, info)')
    .split_multi()
    .min_rep()
    .annotate_variants_expr('va.nAA = (if (va.wasSplit) va.info.Hom_NFE[va.aIndex - 1] else va.info.Hom_NFE[0]).toLong()')
    .annotate_variants_expr('va.nRA = ((if (va.wasSplit) va.info.AC_NFE[va.aIndex - 1] else va.info.AC_NFE[0]) - 2*va.nAA).toLong()')
    .annotate_variants_expr('va.nmissing = (7509 - va.info.AN_NFE/2).toLong()')
    .annotate_variants_expr('va.nRR = (7509 - va.nAA - va.nRA - va.nmissing).toLong()')
    .export_variants('tmp/gnomad_af_table.tsv', 'chrom = v.contig, pos = v.start, ref = v.ref, alt = v.alt(), nRR = va.nRR, nRA = va.nRA, nAA = va.nAA, nmissing = va.nmissing')
)
