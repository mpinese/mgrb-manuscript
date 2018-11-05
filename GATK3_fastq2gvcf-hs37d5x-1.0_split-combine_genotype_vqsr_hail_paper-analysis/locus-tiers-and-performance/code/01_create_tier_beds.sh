#!/bin/bash
set -euo pipefail
set -x

# Rules (from GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqcode/02_variantqc_add_locus_tiers.py):
#
# if xsome ~ "^([0-9]+|X|Y)$": Tier 3
# elif badCoverage_dist <= 5 || badComplexity_dist <= 5 || badMappability_dist <= 5 || badEncodeExcluded_dist <= 5 || va.locus.badPAR: Tier 3
# elif goodGiaBHighConf_dist > 0: Tier 2
# elif xsome ~ "^(X|Y)$": Tier 2
# else: Tier 1
#
# Re-implement on command line with bedtools / awk

mkdir -p tmp

# First define the genome
gzip -dc ../../../locus-annotations/hs37d5x_data/genome.bed.gz > tmp/genome.bed


# TIER 3
######################

# Define canonical xsome tier 3 regions based on an expansion of the bad beds
gzip -dc ../../../locus-annotations/regions/bad.depth.bed.gz | awk 'BEGIN {FS="\t";OFS="\t"} {$2 = $2 - 5; $3 = $3 + 5; if ($2 < 0) $2 = 0; print $1, $2, $3}' | bedtools merge | bedtools intersect -a - -b tmp/genome.bed \
    > tmp/tier3.coverage.extended.bed
gzip -dc ../../../locus-annotations/regions/bad.mdust.bed.gz | awk 'BEGIN {FS="\t";OFS="\t"} {$2 = $2 - 5; $3 = $3 + 5; if ($2 < 0) $2 = 0; print $1, $2, $3}' | bedtools merge | bedtools intersect -a - -b tmp/genome.bed \
    > tmp/tier3.complexity.extended.bed
gzip -dc ../../../locus-annotations/regions/bad.wgEncodeCrgMapabilityAlign100mer.bed.gz | awk 'BEGIN {FS="\t";OFS="\t"} {$2 = $2 - 5; $3 = $3 + 5; if ($2 < 0) $2 = 0; print $1, $2, $3}' | bedtools merge | bedtools intersect -a - -b tmp/genome.bed \
    > tmp/tier3.mappability.extended.bed
gzip -dc ../../../locus-annotations/regions/bad.wgEncodeExcludable.bed.gz | awk 'BEGIN {FS="\t";OFS="\t"} {$2 = $2 - 5; $3 = $3 + 5; if ($2 < 0) $2 = 0; print $1, $2, $3}' | bedtools merge | bedtools intersect -a - -b tmp/genome.bed \
    > tmp/tier3.encodeexcluded.extended.bed

# And additional tier 3 non-canonical contigs
grep -P '^([0-9]+|X|Y)\t' tmp/genome.bed > tmp/genome.canonical.bed
bedtools subtract -a tmp/genome.bed -b tmp/genome.canonical.bed > tmp/tier3.non_canonical.bed

# Combine the above to define tier 3
cat tmp/tier3.coverage.extended.bed tmp/tier3.complexity.extended.bed tmp/tier3.mappability.extended.bed tmp/tier3.encodeexcluded.extended.bed tmp/tier3.non_canonical.bed | sort -k1,1 -k2,2n | bedtools merge > tmp/tier3.bed


# TIER 2
######################

# Sex chromosomes
grep -P '^(X|Y)\t' tmp/genome.bed > tmp/tier23.gonosomes.bed
bedtools subtract -a tmp/tier23.gonosomes.bed -b tmp/tier3.bed > tmp/tier2.gonosomes.bed

# Not in GiaB HC
gzip -dc ../../../locus-annotations/regions/good.HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed.gz | sort -k1,1 -k2,2n | bedtools subtract -a tmp/genome.bed -b - > \
    tmp/tier23.notgiab.bed
bedtools subtract -a tmp/tier23.notgiab.bed -b tmp/tier3.bed > tmp/tier2.notgiab.bed

# Combine to define tier 2
cat tmp/tier2.gonosomes.bed tmp/tier2.notgiab.bed | sort -k1,1 -k2,2n | bedtools merge > tmp/tier2.bed


# TIER 1
#####################
bedtools subtract -a tmp/genome.bed -b tmp/tier3.bed | bedtools subtract -a - -b tmp/tier2.bed > tmp/tier1.bed
