#!./bin/pyhail.sh
import hail

hc = hail.HailContext(log = 'log/04_generate_common_bed.log', tmp_dir = 'tmp/hail')

vds_gnomad_split = hc.read('../gnomad.genomes.r2.0.1.sites.combined.split.minrep.vds')

vds_gnomad_split_common = vds_gnomad_split.filter_variants_expr('''
    v.altAllele().isSNP() &&
    (va.gnomad.AF_AFR >= 0.05 || 
    va.gnomad.AF_AMR >= 0.05 || 
    va.gnomad.AF_ASJ >= 0.05 || 
    va.gnomad.AF_EAS >= 0.05 || 
    va.gnomad.AF_FIN >= 0.05 || 
    va.gnomad.AF_NFE >= 0.05 || 
    va.gnomad.AF_OTH >= 0.05)''')

vds_gnomad_split_common.export_variants('./tmp/04_gnomad.genomes.r2.0.1.sites.combined.split.minrep.common_snps.unsorted_unmerged.bed', 'v.contig, v.start-1, v.start')

# Note: this output file still needs sorting and merging.
