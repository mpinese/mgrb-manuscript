#!/bin/bash

if [ ! -e "../coverage/genomes.done" ]; then
  mkdir -p ../coverage/
  gsutil -m cp -r gs://gnomad-public/release/2.0.1/coverage/genomes/ ../coverage/ && touch ../coverage/genomes.done
fi

function make_coverage_bed
{
  # params: infile depth frac outfile
  gzip -dc "$1" | python3 coverage2bed.py /dev/stdin "$2" "$3" "$4"
}

export -f make_coverage_bed

parallel make_coverage_bed ../coverage/genomes/gnomad.genomes.r2.0.1.chr{}.coverage.txt.gz 15 0.98 ../coverage/gnomad.genomes.r2.0.1.chr{}.15_098.bed ::: $(seq 1 22)

cat ../coverage/gnomad.genomes.r2.0.1.chr*.15_098.bed | sort -k1,1 -k2,2n > ../coverage/gnomad.genomes.r2.0.1.15_098.bed

awk 'BEGIN{FS="\t";OFS="\t"} ($3-$2>5);' < ../coverage/gnomad.genomes.r2.0.1.15_098.bed > ../coverage/gnomad.genomes.r2.0.1.15_098.min5bp.bed

