#!/bin/bash
set -euo pipefail

function convert_ukbb
{
    infreq="$1"
    indbsnp="$2"
    outdb="$3"

    tmpfile=$(mktemp)

    # Note that 3442 loci (of 8034019, < 0.05%) are rejected in the following by the rsid filter.
    xz -dc "${infreq}" | awk 'BEGIN {OFS="\t"} (NR > 1) {print $2, $3, $4, $5, $6, $7, $10}' | grep -P '^rs[0-9]+\t' | sed 's/^rs//' > "${tmpfile}"

    sqlite3 <<EOF
ATTACH "${indbsnp}" AS dbsnp;

CREATE TABLE dbsnp_counts (
  rsid     INTEGER PRIMARY KEY,
  cref     TEXT    NOT NULL,
  calt     TEXT    NOT NULL,
  cRR      INTEGER NOT NULL,
  cRA      INTEGER NOT NULL,
  cAA      INTEGER NOT NULL,
  nmissing INTEGER NOT NULL);
.mode tabs
.import "${tmpfile}" dbsnp_counts

CREATE TABLE merged_counts AS
  SELECT chrom, pos, ref, alt, cref, calt, cRR, cRA, cAA, nmissing
  FROM dbsnp_counts INNER JOIN dbsnp.dbsnp AS dbs ON dbsnp_counts.rsid = dbs.rsid WHERE dbs.multiallelic = 0;

.mode tabs
.header off
.output "${tmpfile}"
SELECT * FROM merged_counts ORDER BY chrom, pos, ref, alt;

EOF

    tmpfile2=$(mktemp)
    awk 'BEGIN {OFS="\t"} ($3 == $5 && $4 == $6) { print $1, $2, $3, $4, $7, $8, $9, $10 } ($3 == $6 && $4 == $5) { print $1, $2, $3, $4, $9, $8, $7, $10 }' < "${tmpfile}" > "${tmpfile2}"

    sqlite3 "${outdb}" <<EOF
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
.import "${tmpfile2}" allelefreqs
EOF

    rm -f "${tmpfile}" "${tmpfile2}"
}


export -f convert_ukbb


parallel convert_ukbb ../UKBB_Sini/GENO_COUNTS/UKB_AFTER_QC_WBritish_GENOCOUNT_AGE_{1}.txt.xz ../../commonvar-burden/dbSNP_All_20180423_strandspecific_biallelic_snvs.db ../UKBB.mf.{1}.db ::: less55 55-60 60-65 65-70 70-75 more75
