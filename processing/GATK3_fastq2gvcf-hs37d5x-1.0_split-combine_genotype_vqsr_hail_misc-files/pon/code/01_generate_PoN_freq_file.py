#!./bin/pyhail.sh
import hail


hc = hail.HailContext(log = 'log/01_generate_PoN_freq_file.log', tmp_dir = 'tmp/hail')

vds_mgrb = hc.read('../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.vds')

# Split variants
vds_mgrb_split = vds_mgrb.split_multi().min_rep()

# Export as-is for Vel
(vds_mgrb_split
    .export_variants(
        '../01_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.PoN-with_singletons.txt', 
        'chrom = v.contig, pos = v.start, ref = v.altAllele().ref, alt = v.altAllele().alt, filters = va.variant.filters.mkString(","), CC = va.variant.metrics.chromatid_counts.total, AC = va.alleles[va.aIndex-1].metrics.allele_counts.total')
)

# Filter to non-singletons and export as PoN file for Marie and Mark C.
#(vds_mgrb_split
#    .filter_variants_expr('va.alleles[va.aIndex-1].metrics.allele_counts.total > 1')
#    .export_variants(
#        '../01_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.PoN-no_singletons.txt', 
#        'chrom = v.contig, pos = v.start, ref = v.altAllele().ref, alt = v.altAllele().alt, filters = va.variant.filters.mkString(","), CC = va.variant.metrics.chromatid_counts.total, AC = va.alleles[va.aIndex-1].metrics.allele_counts.total')
#)
