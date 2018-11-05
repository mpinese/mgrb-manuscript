#!/bin/bash
set -euo pipefail

mkdir -p tmp

gzip -dc ../swegen_20180409/anon-SweGen_STR_NSPHS_1000samples_SNV_hg19.vcf.gz | python3 swegen_splitvcf.py | bgzip > tmp/swegen_split.vcf.bgz

./bin/pyhail.sh 01h_swegen_hail_minrep.py

# SweGen 20180409 seems to have two errors: duplicated variants
# 2:230749380:G:T and 9:140840289:T:TTTA.  In both cases the first
# variant seems to be more reliable (nmissing = 0).
# Explicitly filter out the lower-quality variants.
awk 'BEGIN{FS="\t";OFS="\t";last_var=""} (($1 ":" $2 ":" $3 ":" $4) != last_var) { last_var = $1 ":" $2 ":" $3 ":" $4; print }' < tmp/swegen_split_minrep_autosomes.tsv > tmp/swegen_split_minrep_autosomes_fixed.tsv

sqlite3 ../swegen_20180409.split.minrep.autosomes.allelefreqs.db <<EOF
DROP TABLE IF EXISTS allelefreqs;

CREATE TABLE allelefreqs (
chrom    TEXT    NOT NULL,
pos      INTEGER NOT NULL,
ref      TEXT    NOT NULL,
alt      TEXT    NOT NULL,
nRR      INTEGER NOT NULL,
nRA      INTEGER NOT NULL,
nAA      INTEGER NOT NULL,
nmissing INTEGER NOT NULL,
PRIMARY KEY (chrom, pos, ref, alt));

.mode tabs
.import tmp/swegen_split_minrep_autosomes_fixed.tsv allelefreqs
EOF

rm -f tmp/*

