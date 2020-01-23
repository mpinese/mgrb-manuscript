#!/bin/bash
set -euo pipefail

./bin/pyhail.sh 05h_gnomad2af_tbl.py

function aftbl2db
{
  tbl="$1"
  db="$2"

  tr '\t' ',' < "$tbl" | awk 'BEGIN{FS=",";OFS=","} ($6 != 0 || $7 != 0);' > "${tbl}.csv"

  tail -n+2 "${tbl}.csv" > "${tbl}.noheader.csv"

  sqlite3 "${db}" <<EOF
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
.mode csv
.import "${tbl}.noheader.csv" allelefreqs
EOF
}


aftbl2db ./tmp/gnomad_af_table.tsv ../gnomad.genomes.r2.0.1.sites.autosomes.split.minrep.NFE.allelefreqs.db

# Note: these dbs are missing as GnomAD does not report sex-specific allele counts for subpopulations.
#aftbl2db ./tmp/gnomad_af_table_m.tsv ../gnomad.genomes.r2.0.1.sites.autosomes.split.minrep.NFE.m.allelefreqs.db
#aftbl2db ./tmp/gnomad_af_table_f.tsv ../gnomad.genomes.r2.0.1.sites.autosomes.split.minrep.NFE.f.allelefreqs.db
