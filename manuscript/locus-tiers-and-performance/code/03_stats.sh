#!/bin/bash

echo "Relative to hs37d5x:"
echo "Tier 1:    $(awk 'BEGIN{total=0} {total += $3-$2} END {printf "%10d", total}' < tmp/tier1.bed) bp"
echo "Tier 2:    $(awk 'BEGIN{total=0} {total += $3-$2} END {printf "%10d", total}' < tmp/tier2.bed) bp"
echo "Tier 3:    $(awk 'BEGIN{total=0} {total += $3-$2} END {printf "%10d", total}' < tmp/tier3.bed) bp"
echo "Genome:    $(awk 'BEGIN{total=0} {total += $3-$2} END {printf "%10d", total}' < tmp/genome.bed) bp"
echo
echo "Relative to canonical contigs:"
echo "Tier 1:    $(bedtools intersect -a tmp/tier1.bed -b tmp/genome.canonical.bed | awk 'BEGIN{total=0} {total += $3-$2} END {printf "%10d", total}') bp"
echo "Tier 2:    $(bedtools intersect -a tmp/tier2.bed -b tmp/genome.canonical.bed | awk 'BEGIN{total=0} {total += $3-$2} END {printf "%10d", total}') bp"
echo "Tier 3:    $(bedtools intersect -a tmp/tier3.bed -b tmp/genome.canonical.bed | awk 'BEGIN{total=0} {total += $3-$2} END {printf "%10d", total}') bp"
echo "Canonical: $(awk 'BEGIN{total=0} {total += $3-$2} END {printf "%10d", total}' < tmp/genome.canonical.bed) bp"
echo
echo "Relative to refseq coding:"
echo "Tier 1:    $(bedtools intersect -a tmp/tier1.bed -b tmp/refgene.codingexons.bed | awk 'BEGIN{total=0} {total += $3-$2} END {printf "%10d", total}') bp"
echo "Tier 2:    $(bedtools intersect -a tmp/tier2.bed -b tmp/refgene.codingexons.bed | awk 'BEGIN{total=0} {total += $3-$2} END {printf "%10d", total}') bp"
echo "Tier 3:    $(bedtools intersect -a tmp/tier3.bed -b tmp/refgene.codingexons.bed | awk 'BEGIN{total=0} {total += $3-$2} END {printf "%10d", total}') bp"
echo "Coding:    $(awk 'BEGIN{total=0} {total += $3-$2} END {printf "%10d", total}' < tmp/refgene.codingexons.bed) bp"
