#!/bin/bash
set -euo pipefail

REFERENCE="../../../resources/hs37d5x/reference_genome/hs37d5x.fa"

./bin/mdust "${REFERENCE}" -c | awk 'BEGIN {FS="\t";OFS="\t"} {print $1, $3-1, $4}' | sort -k1,1 -k2,2n | ./bin/bedtools merge | bgzip > ../hs37d5x_data/bad.mdust.bed.gz
