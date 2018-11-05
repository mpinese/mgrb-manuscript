#!/bin/bash
set -euo pipefail

REFERENCE="../../../resources/hs37d5x/reference_genome/hs37d5x.fa"

# Construct a genome-wide .bed.  This includes decoys, patch contigs, etc.
awk 'BEGIN{FS="\t";OFS="\t"} {print $1,0,$2}' < "${REFERENCE}.fai" | sort -k1,1 -k2,2n | bgzip > ../hs37d5x_data/genome.bed.gz

# Construct a bed of canonical chromosomes only.
gzip -dc ../hs37d5x_data/genome.bed.gz | grep -P '^([0-9]+|X|Y|MT)\t' | bgzip > ../hs37d5x_data/genome.canonical.bed.gz

# Define decoy sequences as always bad regions
gzip -dc ../hs37d5x_data/genome.bed.gz | awk 'BEGIN{FS="\t";OFS="\t"} ($1 == "hs37d5" || $1 == "NC_001422.1" || $1 == "NC_007605")' | bgzip > ../hs37d5x_data/bad.decoy.bed.gz

./bin/bigWigToBedGraph ../source_data/wgEncodeCrgMapabilityAlign100mer.bigWig /dev/stdout | sed 's/^chr//; s/^M/MT/' | awk 'BEGIN{FS="\t";OFS="\t"} ($4 < 1) {print $1,$2,$3}' | sort -k1,1 -k2,2n | ./bin/bedtools merge | bgzip > ../hs37d5x_data/bad.wgEncodeCrgMapabilityAlign100mer.bed.gz
./bin/bigWigToBedGraph ../source_data/wgEncodeCrgMapabilityAlign100mer.bigWig /dev/stdout | sed 's/^chr//; s/^M/MT/' | awk 'BEGIN{FS="\t";OFS="\t"} ($4 == 1) {print $1,$2,$3}' | sort -k1,1 -k2,2n | ./bin/bedtools merge | bgzip > ../hs37d5x_data/good.wgEncodeCrgMapabilityAlign100mer.bed.gz

./bin/bigWigToBedGraph ../source_data/wgEncodeDukeMapabilityUniqueness35bp.bigWig /dev/stdout | sed 's/^chr//; s/^M/MT/' | awk 'BEGIN{FS="\t";OFS="\t"} ($4 < 1) {print $1,$2,$3}' | sort -k1,1 -k2,2n | ./bin/bedtools merge | bgzip > ../hs37d5x_data/bad.wgEncodeDukeMapabilityUniqueness35bp.bed.gz
./bin/bigWigToBedGraph ../source_data/wgEncodeDukeMapabilityUniqueness35bp.bigWig /dev/stdout | sed 's/^chr//; s/^M/MT/' | awk 'BEGIN{FS="\t";OFS="\t"} ($4 == 1) {print $1,$2,$3}' | sort -k1,1 -k2,2n | ./bin/bedtools merge | bgzip > ../hs37d5x_data/good.wgEncodeDukeMapabilityUniqueness35bp.bed.gz

gzip -dc ../source_data/wgEncodeDacMapabilityConsensusExcludable.bed.gz | sed 's/^chr//; s/^M/MT/' | sort -k1,1 -k2,2n | ./bin/bedtools merge | bgzip > ../hs37d5x_data/bad.wgEncodeDacMapabilityConsensusExcludable.bed.gz

gzip -dc ../source_data/wgEncodeDukeMapabilityRegionsExcludable.bed.gz | sed 's/^chr//; s/^M/MT/' | sort -k1,1 -k2,2n | ./bin/bedtools merge | bgzip > ../hs37d5x_data/bad.wgEncodeDukeMapabilityRegionsExcludable.bed.gz

gzip -dc ../source_data/rmsk.txt.gz | cut -f 6-8 | sed 's/^chr//; s/^M/MT/' | sort -k1,1 -k2,2n | ./bin/bedtools merge | bgzip > ../hs37d5x_data/bad.rmsk.bed.gz

sort -k1,1 -k2,2n ../source_data/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed | bgzip > ../hs37d5x_data/good.HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed.gz
