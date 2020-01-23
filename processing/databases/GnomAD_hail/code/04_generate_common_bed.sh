#!/bin/bash
set -euo pipefail

./bin/pyhail.sh 04h_generate_common_bed.py

sort -k1,1 -k2,2n tmp/04_gnomad.genomes.r2.0.1.sites.combined.split.minrep.common_snps.unsorted_unmerged.bed | ./bin/bedtools merge | bgzip > ../04_gnomad.genomes.r2.0.1.sites.combined.split.minrep.common_snps.bed.gz
