#!/bin/bash
set -euo pipefail

mkdir -p temp


export SQLITE_TMPDIR=./temp


# Extract basic information from dbSNP 151 (human_9606_b151_GRCh37p13) for subsequent import to SQL
gzip -dc dbsnp_All_20180423.vcf.gz \
  | sed 's/^chr//' \
  | awk '
BEGIN {
  FS="\t";OFS=","
}

/^[^#]/ {
  i++
  chrom = $1
  pos = $2
  ref = $4
  alts = $5
  info = $8

  nalts = split(alts, alts_split, ",")

  # Drop multi-allelics, indels/MNVs, and non strand-specific SNPs
  # Also drop anything not autosomal
  # Strand specific enforced at the individual alt level below
  multiallelic = nalts > 1
  notsnv = length(ref) != 1 || length(alts) != 1
  autosomal = match(chrom, "^[0-9]+$") != 0

  if (multiallelic || notsnv || !autosomal)
    next

  ninfo = split(info, info_split, ";")
  negstrand = 0
  rsid = "ERROR"
  for (info_field_i in info_split) {
    info_field = info_split[info_field_i]
    info_field_size = split(info_field, info_field_split, "=")
    if (info_field_size == 1) {
      if (info_field == "RV")
        negstrand = 1
    } else {
      if (info_field_split[1] == "RS")
        rsid = info_field_split[2]
    }
  }

  j++
  k += nalts
  for (alt_i in alts_split) {
    alt = alts_split[alt_i]
    strand_specific = ((ref == "G" || ref == "C") && (alt == "A" || alt == "T")) || ((ref == "A" || ref == "T") && (alt == "G" || alt == "C"))
    if (strand_specific) {
      print rsid, chrom, pos, ref, alt, multiallelic, negstrand
    }
  }

  if (i >= nextprint_i) {
    print chrom, i, j, k > "/dev/stderr"
    nextprint_i = i + 1000000
  }
}' \
  > temp/dbsnp_sssnv.csv


rm -f ../dbSNP_All_20180423_strandspecific_biallelic_snvs.db
sqlite3 ../dbSNP_All_20180423_strandspecific_biallelic_snvs.db <<EOF
CREATE TABLE dbsnp (
    rsid         INTEGER PRIMARY KEY,
    chrom        TEXT,
    pos          INTEGER,
    ref          TEXT,
    alt          TEXT,
    multiallelic INTEGER,
    negstrand    INTEGER
);
.mode csv
.import temp/dbsnp_sssnv.csv dbsnp
CREATE INDEX dbsnp_idx_vid ON dbsnp(chrom, pos, ref, alt);

EOF


# Extract the rsids for the EBI GWAS loci
cut -f 24 gwas_catalogue.tsv | tail -n +2 | sort | uniq | grep -E '^[0-9]+$' > temp/ebi_snpids.csv

# Extract chrom, pos, ref, alt for the manually-entered models
# cut -f2 manual_polygenic_scores.GnomAD_NFE_AFs.models | tail -n+2 | grep -v '^OFFSET$' | sort | uniq | sed 's/:/,/g' > temp/grs_variants.csv
cut -f2 -d, manual_polygenic_scores.hcr_tag_rescued.csv | tail -n+2 | grep -v '^OFFSET$' | sort | uniq | sed 's/:/,/g' > temp/grs_variants.csv

# Import dbSNP and the EBI GWAS loci to SQL
# Intersect with the AFs

sqlite3 <<EOF

ATTACH "../../allelefreqs/MGRB_GnomAD_SweGen_cancer_UKBB_AFs_outerjoined.db" AS aforig;
ATTACH "../dbSNP_All_20180423_strandspecific_biallelic_snvs.db" AS dbsnp;

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
.output ../MGRB_GnomAD_SweGen_cancer_UKBB_AFs_outerjoined_ss-auto-gwas-snps.csv
SELECT * FROM ebi_grs_dbsnp_afwide WHERE chrom IS NOT NULL ORDER BY chrom ASC, pos ASC, ref ASC, alt ASC;

.output ../MGRB_GnomAD_SweGen_cancer_UKBB_AFs_outerjoined_ss-auto-gwas-snps_hcr.csv
SELECT * FROM ebi_grs_dbsnp_afwide_highconf WHERE chrom IS NOT NULL ORDER BY chrom ASC, pos ASC, ref ASC, alt ASC;

EOF


rm -f ../MGRB_GnomAD_SweGen_cancer_UKBB_AFs_outerjoined_ss-auto-gwas-snps.csv.xz
xz -9 ../MGRB_GnomAD_SweGen_cancer_UKBB_AFs_outerjoined_ss-auto-gwas-snps.csv

rm -f ../MGRB_GnomAD_SweGen_cancer_UKBB_AFs_outerjoined_ss-auto-gwas-snps_hcr.csv.xz
xz -9 ../MGRB_GnomAD_SweGen_cancer_UKBB_AFs_outerjoined_ss-auto-gwas-snps_hcr.csv
