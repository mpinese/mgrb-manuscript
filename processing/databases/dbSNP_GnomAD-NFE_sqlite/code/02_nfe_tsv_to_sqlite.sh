#!/bin/bash
set -euo pipefail

sqlite3 ../02_gnomad.genomes.r2.0.1.sites.combined.split.minrep.NFE.dbSNP.autosomes.dba <<END_OF_SCRIPT
DROP TABLE IF EXISTS dbsnp;
CREATE TABLE dbsnp (
  chrom  TEXT     NOT NULL,
  pos    INTEGER  NOT NULL,
  ref    TEXT     NOT NULL,
  alt    TEXT     NOT NULL,
  rsid   TEXT     NOT NULL,
  AC     INTEGER  NOT NULL,
  AN     INTEGER  NOT NULL,
  PRIMARY KEY (chrom, pos, ref, alt));
CREATE INDEX dbsnp_rsid_idx ON dbsnp(rsid);
.separator "\t"
.import ../01_gnomad.genomes.r2.0.1.sites.combined.split.minrep.NFE.dbSNP.autosomes.tsv dbsnp
END_OF_SCRIPT
