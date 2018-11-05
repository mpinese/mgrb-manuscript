#!/bin/bash
set -euo pipefail

mkdir -p ../regions

# Tier QC (used for QC and comparative work):
# (Autosomes intersect GiaB gold standard regions)
# - mdust
# - rmsk
# - wgEncodeCrgMapabilityAlign100mer score < 1
# - wgEncodeDacMapabilityConsensusExcludable
# - wgEncodeDukeMapabilityRegionsExcludable
# - wgEncodeDukeMapabilityUniqueness35bp score < 1
#
# Interpretation: this tier is based on various
# sequence-only and database resources (ie not built
# from observed properties of the MGRB data).  It is
# a very conservative filter intended to capture very
# high confidence regions for QC.  Good data should
# have excellent performance in these regions.
#
# Use a morphological opening operation at the end
# to remove runs of < 100 bp.  Rough testing indicates
# that this leads to the loss of ~6% of bases in regions.
#
# Follow up with a final erosion by 10 bp, to ensure
# regions do not sit on the edge of good sequence.
#
# The resulting qc.bed.gz contains 1131171577 bp,
# about 40% of the total autosome size of 2881033286 bp.

./bin/bedtools intersect -a ../hs37d5x_data/genome.canonical.bed.gz -b ../hs37d5x_data/good.HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed.gz | \
  grep -P '^[0-9]+\t' | \
  ./bin/bedtools subtract -a - -b ../hs37d5x_data/bad.mdust.bed.gz | \
  ./bin/bedtools subtract -a - -b ../hs37d5x_data/bad.rmsk.bed.gz | \
  ./bin/bedtools subtract -a - -b ../hs37d5x_data/bad.wgEncodeCrgMapabilityAlign100mer.bed.gz | \
  ./bin/bedtools subtract -a - -b ../hs37d5x_data/bad.wgEncodeDacMapabilityConsensusExcludable.bed.gz | \
  ./bin/bedtools subtract -a - -b ../hs37d5x_data/bad.wgEncodeDukeMapabilityRegionsExcludable.bed.gz | \
  ./bin/bedtools subtract -a - -b ../hs37d5x_data/bad.wgEncodeDukeMapabilityUniqueness35bp.bed.gz | \
  awk 'BEGIN{FS="\t";OFS="\t"} ($3-$2 >= 100);' | \
  awk 'BEGIN{FS="\t";OFS="\t"} {$2 += 10; $3 -= 10; print}' | \
  sort -k1,1 -k2,2n | \
  bgzip > ../regions/qc.bed.gz


# Tier LoH (used for subclonal LoH detection).
# As for Tier QC, plus poor depth regions removed,
# finished with a final opening of 100 bp and
# erosion of 10 bp.
./bin/bedtools intersect -a ../hs37d5x_data/genome.canonical.bed.gz -b ../hs37d5x_data/good.HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed.gz | \
  grep -P '^[0-9]+\t' | \
  ./bin/bedtools subtract -a - -b ../hs37d5x_data/bad.mdust.bed.gz | \
  ./bin/bedtools subtract -a - -b ../hs37d5x_data/bad.rmsk.bed.gz | \
  ./bin/bedtools subtract -a - -b ../hs37d5x_data/bad.wgEncodeCrgMapabilityAlign100mer.bed.gz | \
  ./bin/bedtools subtract -a - -b ../hs37d5x_data/bad.wgEncodeDacMapabilityConsensusExcludable.bed.gz | \
  ./bin/bedtools subtract -a - -b ../hs37d5x_data/bad.wgEncodeDukeMapabilityRegionsExcludable.bed.gz | \
  ./bin/bedtools subtract -a - -b ../hs37d5x_data/bad.wgEncodeDukeMapabilityUniqueness35bp.bed.gz | \
  ./bin/bedtools subtract -a - -b <(gzip -dc ../gooddepth_samplecounts.bed.gz | awk "BEGIN{FS=\"\t\";OFS=\"\t\"} (\$4 / $(cat ../gooddepth_totalsamples.txt) < 0.5);") | \
  awk 'BEGIN{FS="\t";OFS="\t"} ($3-$2 >= 100);' | \
  awk 'BEGIN{FS="\t";OFS="\t"} {$2 += 10; $3 -= 10; print}' | \
  sort -k1,1 -k2,2n | \
  bgzip > ../regions/loh.bed.gz
bedtools subtract -a ../hs37d5x_data/genome.bed.gz -b ../regions/loh.bed.gz | bgzip > ../regions/loh-blacklist.bed.gz

# Tier CNV (used for CN comparisons, no GC exclusion,
# no moderately poor mappability exclusion)
# Assumes that the downstream CN tools can perform
# GC and simple mappability correction.
#
# Canonical genome
# - mdust
# - rmsk
# - wgEncodeCrgMapabilityAlign100mer score < 1
# - wgEncodeDacMapabilityConsensusExcludable
# - wgEncodeDukeMapabilityRegionsExcludable
# - regions with unusual empirical depth in >= 50% of the cohort.
# dilation, kernel size 10
# erosion, kernel size 450
./bin/bedtools subtract -a ../hs37d5x_data/genome.canonical.bed.gz -b ../hs37d5x_data/bad.mdust.bed.gz | \
  ./bin/bedtools subtract -a - -b ../hs37d5x_data/bad.rmsk.bed.gz | \
  ./bin/bedtools subtract -a - -b ../hs37d5x_data/bad.wgEncodeCrgMapabilityAlign100mer.bed.gz | \
  ./bin/bedtools subtract -a - -b ../hs37d5x_data/bad.wgEncodeDacMapabilityConsensusExcludable.bed.gz | \
  ./bin/bedtools subtract -a - -b ../hs37d5x_data/bad.wgEncodeDukeMapabilityRegionsExcludable.bed.gz | \
  ./bin/bedtools subtract -a - -b <(gzip -dc ../gooddepth_samplecounts.bed.gz | awk "BEGIN{FS=\"\t\";OFS=\"\t\"} (\$4 / $(cat ../gooddepth_totalsamples.txt) < 0.5);")
  awk 'BEGIN{FS="\t";OFS="\t"} {$2 -= 10/2; $2 += 10/2; if ($2<0) {$2=0}; print}' | \
  ./bin/bedtools merge | \
  awk 'BEGIN{FS="\t";OFS="\t"} ($3 - $2 > 450);' | \
  bgzip > ../regions/cnv.bed.gz


# Create a combined ENCODE excludable track
cat ../hs37d5x_data/bad.wgEncodeDacMapabilityConsensusExcludable.bed.gz ../hs37d5x_data/bad.wgEncodeDukeMapabilityRegionsExcludable.bed.gz | \
  gzip -dc | \
  sort -k1,1 -k2,2n | \
  ./bin/bedtools merge | \
  ./bin/bedtools intersect -a - -b ../hs37d5x_data/genome.bed.gz | \
  bgzip > ../regions/bad.wgEncodeExcludable.bed.gz


# Intersect the UCSC-derived tracks with the genome, to remove ucsc-specific contigs
./bin/bedtools intersect -a ../hs37d5x_data/bad.rmsk.bed.gz -b ../hs37d5x_data/genome.bed.gz | bgzip > ../regions/bad.rmsk.bed.gz
./bin/bedtools intersect -a ../hs37d5x_data/bad.wgEncodeCrgMapabilityAlign100mer.bed.gz -b ../hs37d5x_data/genome.bed.gz | bgzip > ../regions/bad.wgEncodeCrgMapabilityAlign100mer.bed.gz


# Copy over other regions unchanged
cp ../hs37d5x_data/bad.depth.bed.gz ../regions/
cp ../hs37d5x_data/bad.mdust.bed.gz ../regions/
cp ../hs37d5x_data/good.HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed.gz ../regions/

# Index everything
for f in ../regions/*.bed.gz; do tabix $f; done
