#!/bin/bash
set -euo pipefail

# GnomAD AFs are in the database
# ../gnomad.genomes.r2.0.1.sites.combined.split.minrep.NFE.dbSNP.autosomes.dba
# with schema:
# CREATE TABLE dbsnp (
#   chrom  TEXT     NOT NULL,
#   pos    INTEGER  NOT NULL,
#   ref    TEXT     NOT NULL,
#   alt    TEXT     NOT NULL,
#   rsid   TEXT     NOT NULL,
#   AC     INTEGER  NOT NULL,
#   AN     INTEGER  NOT NULL,
#   PRIMARY KEY (chrom, pos, ref, alt));
# CREATE INDEX dbsnp_rsid_idx ON dbsnp(rsid);

# MGRB AFs are in the database
# ../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.popstar.allele_freqs.dba
# with schema:
# CREATE TABLE mgrb (
#   chrom    TEXT     NOT NULL,
#   pos      INTEGER  NOT NULL,
#   ref      TEXT     NOT NULL,
#   alt      TEXT     NOT NULL,
#   AC_MGRB  INTEGER  NOT NULL,
#   AN_MGRB  INTEGER  NOT NULL,
#   PRIMARY KEY (chrom, pos, ref, alt));


# Perform an outer join of these two tables to form a combined AF db.
sqlite3 ../MGRB_GnomAD_AFs_outerjoined.dba <<EOF
ATTACH "../gnomad.genomes.r2.0.1.sites.combined.split.minrep.NFE.dbSNP.autosomes.dba" AS gnomad;
ATTACH "../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.popstar.allele_freqs.dba" AS mgrb;

CREATE TABLE mgrb_gnomad AS
  SELECT g.chrom, g.pos, g.ref, g.alt, g.AC AS AC_gnomad, g.AN AS AN_gnomad, m.AC_mgrb, m.AN_mgrb
    FROM gnomad.dbsnp g LEFT JOIN mgrb.mgrb m
    ON g.chrom = m.chrom AND g.pos = m.pos AND g.ref = m.ref and g.alt = m.alt
  UNION ALL
  SELECT m.chrom, m.pos, m.ref, m.alt, g.AC AS AC_gnomad, g.AN AS AN_gnomad, m.AC_MGRB, m.AN_MGRB
    FROM mgrb.mgrb m LEFT JOIN gnomad.dbsnp g
    ON g.chrom = m.chrom AND g.pos = m.pos AND g.ref = m.ref and g.alt = m.alt
    WHERE g.chrom IS NULL;

DETACH mgrb;
DETACH gnomad;
EOF


# Export common variants in both cohorts as a CSV
sqlite3 ../MGRB_GnomAD_AFs_outerjoined.dba <<EOF
.headers on
.mode csv
.output ../MGRB_GnomAD_AFs_outerjoined_common.csv
SELECT * FROM mgrb_gnomad WHERE AN_gnomad NOT NULL AND AN_mgrb NOT NULL AND 1000*AC_mgrb > AN_mgrb AND 1000*AC_gnomad > AN_gnomad;
EOF

