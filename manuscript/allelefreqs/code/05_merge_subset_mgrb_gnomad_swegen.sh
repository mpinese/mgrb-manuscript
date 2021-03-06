#!/bin/bash
set -euo pipefail
set -x

# MGRB dbs:
# ../MGRB.phase2final.GiaB_HCR.split.minrep.forrvas.dba
# ../MGRB.phase2final.GiaB_HCR.split.minrep.somaticfiltered.forrvas.dba

# GnomAD db:
# ../../../databases/GnomAD_hail/gnomad.genomes.r2.0.1.sites.autosomes.split.minrep.NFE.allelefreqs.dba

# SweGen db:
# ../swegen_20180409.split.minrep.autosomes.allelefreqs.db

# Common schema:
# CREATE TABLE allelefreqs (
#    chrom    TEXT    NOT NULL,
#    pos      INTEGER NOT NULL,
#    ref      TEXT    NOT NULL,
#    alt      TEXT    NOT NULL,
#    nRR      INTEGER NOT NULL,
#    nRA      INTEGER NOT NULL,
#    nAA      INTEGER NOT NULL,
#    nmissing INTEGER NOT NULL,
#    PRIMARY KEY (chrom, pos, ref, alt));

# Perform an outer join of all tables to form a combined AF db.
# This sqlite code is generated by 05h_generate_join_query.py
rm -f ../MGRB_GnomAD_SweGen_cancer_UKBB_AFs_outerjoined.db


mkdir -p temp
export SQLITE_TMPDIR=temp

python3 05h_generate_join_query.py | sqlite3 ../MGRB_GnomAD_SweGen_cancer_UKBB_AFs_outerjoined.db


# Export as a CSV for subsequent bedtools intersection
mkdir -p temp
sqlite3 ../MGRB_GnomAD_SweGen_cancer_UKBB_AFs_outerjoined.db <<EOF
.headers off
.mode csv
.output temp/MGRB_GnomAD_SweGen_cancer_UKBB_AFs_outerjoined.csv
SELECT chrom, pos, ref, alt FROM afwide ORDER BY chrom, pos, alt;
.quit
EOF


# Convert the CSV to a BED
awk 'BEGIN{FS=",";OFS="\t"} {print $1, $2-1, $2+length($3)-1, $0}' < temp/MGRB_GnomAD_SweGen_cancer_UKBB_AFs_outerjoined.csv \
  > temp/MGRB_GnomAD_SweGen_cancer_UKBB_AFs_outerjoined.bed

# Verify correct sort order of the bed
cut -f1 temp/MGRB_GnomAD_SweGen_cancer_UKBB_AFs_outerjoined.bed | uniq > temp/chroms.txt
sort -k1,1 temp/chroms.txt > temp/chroms_sorted.txt
if cmp --silent temp/chroms.txt temp/chroms_sorted.txt; then
  echo "BED chromosome order OK"
else
  exit "BED chromosome order incorrect"
fi

# bedtools intersect with the analysis regions
# Select only variants that are entirely within analysis regions (-f 1)
bedtools intersect -wa -f 1 -a temp/MGRB_GnomAD_SweGen_cancer_UKBB_AFs_outerjoined.bed -b ../MGRBphase2final_dpge15_98pct.GnomAD201_dpge15_98pct.GiaB332HC.bed \
  | cut -f4 \
  > temp/MGRB_GnomAD_SweGen_cancer_UKBB_AFs_sharedhcregions_variants.csv
bedtools intersect -wa -f 1 -a temp/MGRB_GnomAD_SweGen_cancer_UKBB_AFs_outerjoined.bed -b ../MGRBphase2final_dpge15_98pct.GnomAD201_dpge15_98pct.GiaB332HC.ccdsexonsplus2bp.bed \
  | cut -f4 \
  > temp/MGRB_GnomAD_SweGen_cancer_UKBB_AFs_sharedhcregions_ccdsexonsplus2bp_variants.csv


# Re-import these variants back to the master AF db.  They will be used
# to select loci from the main afwide table that are in the respective regions.
sqlite3 ../MGRB_GnomAD_SweGen_cancer_UKBB_AFs_outerjoined.db <<EOF

CREATE TABLE loci_highconf (
  chrom  TEXT,
  pos    INTEGER,
  ref    TEXT,
  alt    TEXT);
.mode csv
.import temp/MGRB_GnomAD_SweGen_cancer_UKBB_AFs_sharedhcregions_variants.csv loci_highconf
CREATE UNIQUE INDEX loci_highconf_idx ON loci_highconf(chrom, pos, ref, alt);

CREATE TABLE loci_highconfcoding (
  chrom  TEXT,
  pos    INTEGER,
  ref    TEXT,
  alt    TEXT);
.mode csv
.import temp/MGRB_GnomAD_SweGen_cancer_UKBB_AFs_sharedhcregions_ccdsexonsplus2bp_variants.csv loci_highconfcoding
CREATE UNIQUE INDEX loci_highconfcoding_idx ON loci_highconfcoding(chrom, pos, ref, alt);

EOF
