#!/bin/bash
set -euo pipefail

mkdir -p ../source_data
cd ../source_data
wget -c http://krishna.gs.washington.edu/download/CADD/v1.3/1000G_phase3.tsv.gz
wget -c http://krishna.gs.washington.edu/download/CADD/v1.3/ESP6500SI.tsv.gz
wget -c http://krishna.gs.washington.edu/download/CADD/v1.3/ExAC_r0.3.tsv.gz
wget -c http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs.tsv.gz
wget -c http://krishna.gs.washington.edu/download/CADD/v1.3/1000G_phase3.tsv.gz.tbi
wget -c http://krishna.gs.washington.edu/download/CADD/v1.3/ESP6500SI.tsv.gz.tbi
wget -c http://krishna.gs.washington.edu/download/CADD/v1.3/ExAC_r0.3.tsv.gz.tbi
wget -c http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs.tsv.gz.tbi
cd -

