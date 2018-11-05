#!/bin/bash

# Comparison regions are:
# * GnomAD >= 98% of samples with DP >= 15, regions of > 5 bp only
#   AND
# * MGRB phase 2 >= 98% of samples with DP >= 15, regions of > 5 bp only
#   AND
# * GiaB gold-standard
#   AND
# * Coding exon +/- 2 bp

# GnomAD regions are at:
# ../../../databases/GnomAD_hail/coverage/gnomad.genomes.r2.0.1.15_098.min5bp.bed

# GiaB gold-standard is at:
# ../../../locus-annotations/source_data/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed

# Coding regions are from UCSC ccdsGene table, exons plus 2 bp, downloaded 20180625 AEST:
gzip -dc ucsc_ccdsGene_exonsplus2bp.bed.gz | cut -f 1-3 | sed 's/^chr//' | sort -k1,1 -k2,2n | bedtools merge > ../ucsc_ccdsGene_exonsplus2bp.hg37.bed

# MGRB regions are generated here:
mkdir -p temp
comm -23 \
    <(cat ../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/tier1.sample_list ../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/tier2.sample_list | grep '^[AB]' | sort) \
    <(cat \
      ../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/45andup_followup_cancer.sample_list \
      ../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/45andup_followup_under70.sample_list \
      ../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.ldpruned.pcrelate.samples_to_drop | sort) \
  > temp/mgrb_phase2_samples.txt

grep -P '^([0-9]+|X|Y)\t' ../../../../resources/hs37d5x/reference_genome/hs37d5x.fa.fai | cut -f 1,2 > temp/hs37d5x.genome

mkdir -p temp/sort
while read sampleid; do
  echo "$sampleid" > /dev/stderr
  gzip -dc ../../../GATK3_fastq2gvcf-hs37d5x-1.0_depth_stats/${sampleid}.dpgood.bed.gz
done < temp/mgrb_phase2_samples.txt \
  | grep -P '^[0-9]+\t' \
  | LC_ALL=C sort -k1,1 -k2,2n -T temp/sort --compress-program=lz4 \
  | bedtools genomecov -g temp/hs37d5x.genome -bg -i /dev/stdin \
  > temp/mgrb_phase2_genomecov.bed
rm -f temp/sort/*

# >= 2519 is >= 98% of MGRB's 2570 samples
awk 'BEGIN{FS="\t";OFS="\t"} ($4 >= 2519);' < temp/mgrb_phase2_genomecov.bed \
  | cut -f1-3 \
  | bedtools merge \
  | awk 'BEGIN{FS="\t";OFS="\t"} ($3-$2 > 5);' \
  > ../MGRBphase2final.dpge15.98pct.min5bp.bed


# Note that MGRB tier 1 was not used as it might be expected to bias
# variant calls away from MGRB, if the equivalent "tier 1" file for GnomAD
# contains different sites.

# The one major uncontrolled variable is whether the GnomAD VDS which was used
# as the ultimate variant input source has bespoke filtering in place.  If so,
# this would be expected to bias variant calls away from GnomAD, for a similar
# reason to the MGRB bias with tier 1 above.

# Intersect
bedtools intersect -a ../../../databases/GnomAD_hail/coverage/gnomad.genomes.r2.0.1.15_098.min5bp.bed -b ../MGRBphase2final.dpge15.98pct.min5bp.bed | \
  bedtools intersect -a - -b ../../../locus-annotations/source_data/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed \
  > ../MGRBphase2final_dpge15_98pct.GnomAD201_dpge15_98pct.GiaB332HC.bed

bedtools intersect -a ../../../locus-annotations/source_data/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed -b ../ucsc_ccdsGene_exonsplus2bp.hg37.bed \
  > ../MGRBphase2final_dpge15_98pct.GnomAD201_dpge15_98pct.GiaB332HC.ccdsexonsplus2bp.bed

