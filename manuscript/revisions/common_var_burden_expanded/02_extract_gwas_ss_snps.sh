#!/bin/bash
set -euo pipefail

mkdir -p temp

export SQLITE_TMPDIR=./temp


# Extract the rsids for the EBI GWAS loci
cut -f 24 ../data/gwas_catalogue.tsv | tail -n +2 | sort | uniq | grep -E '^[0-9]+$' > temp/ebi_snpids.csv

# Extract chrom, pos, ref, alt for the manually-entered models
# cut -f2 manual_polygenic_scores.GnomAD_NFE_AFs.models | tail -n+2 | grep -v '^OFFSET$' | sort | uniq | sed 's/:/,/g' > temp/grs_variants.csv
cut -f2 -d, ../data/manual_polygenic_scores.hcr_tag_rescued.csv | tail -n+2 | grep -v '^OFFSET$' | sort | uniq | sed 's/:/,/g' > temp/grs_variants.csv

# Import dbSNP and the EBI GWAS loci to SQL
# Intersect with the AFs

sqlite3 <<EOF

ATTACH "../data/MGRB_GnomAD_SweGen_cancer_UKBB_AFs_outerjoined.db" AS aforig;
ATTACH "../data/dbSNP_All_20180423_strandspecific_biallelic_snvs.db" AS dbsnp;

-- Load dbSNP IDs for EBI GWAS database variants
CREATE TABLE ebi (rsid INTEGER PRIMARY KEY);
.mode csv
.import temp/ebi_snpids.csv ebi

-- Load vids for manually-entered GRS variants
CREATE TABLE grs (
    chrom        TEXT,
    pos          INTEGER,
    ref          TEXT,
    alt          TEXT
);
.mode csv
.import temp/grs_variants.csv grs
CREATE UNIQUE INDEX grs_idx ON grs(chrom, pos, ref, alt);

-- Extract EBI GWAS database dbSNP info
CREATE TABLE ebi_dbsnp AS
SELECT dbsnp.dbsnp.rsid as rsid, chrom, pos, ref, alt, negstrand
FROM dbsnp.dbsnp INNER JOIN ebi ON dbsnp.dbsnp.rsid = ebi.rsid;

-- Extract manual GRS dbSNP info
CREATE TABLE grs_dbsnp AS
SELECT rsid, g.chrom, g.pos, g.ref, g.alt, negstrand
FROM grs g LEFT JOIN dbsnp.dbsnp d ON g.chrom = d.chrom AND g.pos = d.pos AND g.ref = d.ref AND g.alt = d.alt;

-- Merge the EBI GWAS and manually-entered GRS info
CREATE TABLE ebi_grs_dbsnp AS SELECT * FROM ebi_dbsnp UNION SELECT * FROM grs_dbsnp;
CREATE INDEX ebi_grs_dbsnp_idx ON ebi_grs_dbsnp(chrom, pos, ref, alt);

-- Extract allele freqs for the EBI + GRS loci
CREATE TABLE ebi_grs_dbsnp_afwide AS
SELECT rsid, negstrand, afw.*
FROM ebi_grs_dbsnp vars LEFT JOIN aforig.afwide afw ON
  vars.chrom = afw.chrom AND
  vars.pos = afw.pos AND
  vars.ref = afw.ref AND
  vars.alt = afw.alt;
CREATE UNIQUE INDEX ebi_grs_dbsnp_afwide_idx_coord ON ebi_grs_dbsnp_afwide(chrom, pos, ref, alt);
--CREATE UNIQUE INDEX ebi_grs_dbsnp_afwide_idx_rsid ON ebi_grs_dbsnp_afwide(rsid);

-- Subset the AFs to only high-confidence call regions
CREATE TABLE ebi_grs_dbsnp_afwide_highconf AS
SELECT afw.*
FROM ebi_grs_dbsnp_afwide afw INNER JOIN aforig.loci_highconf hc ON
  afw.chrom = hc.chrom AND
  afw.pos = hc.pos AND
  afw.ref = hc.ref AND
  afw.alt = hc.alt;

.headers on
.mode csv
.output MGRB_GnomAD_SweGen_cancer_UKBB_AFs_outerjoined_ss-auto-gwas-snps.csv
SELECT * FROM ebi_grs_dbsnp_afwide WHERE chrom IS NOT NULL ORDER BY chrom ASC, pos ASC, ref ASC, alt ASC;

.output MGRB_GnomAD_SweGen_cancer_UKBB_AFs_outerjoined_ss-auto-gwas-snps_hcr.csv
SELECT * FROM ebi_grs_dbsnp_afwide_highconf WHERE chrom IS NOT NULL ORDER BY chrom ASC, pos ASC, ref ASC, alt ASC;

EOF


rm -f MGRB_GnomAD_SweGen_cancer_UKBB_AFs_outerjoined_ss-auto-gwas-snps.csv.xz
xz -9 MGRB_GnomAD_SweGen_cancer_UKBB_AFs_outerjoined_ss-auto-gwas-snps.csv

rm -f MGRB_GnomAD_SweGen_cancer_UKBB_AFs_outerjoined_ss-auto-gwas-snps_hcr.csv.xz
xz -9 MGRB_GnomAD_SweGen_cancer_UKBB_AFs_outerjoined_ss-auto-gwas-snps_hcr.csv
