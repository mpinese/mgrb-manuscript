#!/bin/bash

mysql -NB -h genome-mysql.cse.ucsc.edu hg19 -u genome -A -e "SELECT chrom, cdsStart, cdsEnd FROM refGene;" | sed 's/^chr//; s/^M/MT/' | sort -k1,1 -k2,2n | bedtools intersect -a - -b tmp/genome.bed > tmp/refgene.cds.bed
mysql -NB -h genome-mysql.cse.ucsc.edu hg19 -u genome -A -e "SELECT chrom, exonStarts, exonEnds FROM refGene;" | \
    awk 'BEGIN {FS="\t";OFS="\t"} { n=split($2,starts,","); split($3,ends,","); for (i=1;i<n;i++) { print $1, starts[i], ends[i]}}' | \
    sed 's/^chr//; s/^M/MT/' | sort -k1,1 -k2,2n | bedtools intersect -a - -b tmp/genome.bed \
    > tmp/refgene.exons.bed

bedtools intersect -a tmp/refgene.cds.bed -b tmp/refgene.exons.bed > tmp/refgene.codingexons.bed
