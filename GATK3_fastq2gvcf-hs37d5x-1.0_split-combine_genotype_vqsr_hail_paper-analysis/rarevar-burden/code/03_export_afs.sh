#!/bin/bash
set -euo pipefail

sqlite3 <<EOF
ATTACH "../../allelefreqs/MGRB_GnomAD_SweGen_cancer_AFs_outerjoined.db" AS afwide;

.mode csv
.headers on
.output ../MGRB_GnomAD_SweGen_cancer_AFs.MGRBphase2final_dpge15_98pct.GnomAD201_dpge15_98pct.GiaB332HC.ccdsexonsplus2bp.csv
SELECT * FROM afwide.afwide af INNER JOIN afwide.loci_highconfcoding loci
ON af.chrom = loci.chrom AND af.pos = loci.pos AND af.ref = loci.ref AND af.alt = loci.alt;
EOF

rm -f ../MGRB_GnomAD_SweGen_cancer_AFs.MGRBphase2final_dpge15_98pct.GnomAD201_dpge15_98pct.GiaB332HC.ccdsexonsplus2bp.csv.xz
xz -9 ../MGRB_GnomAD_SweGen_cancer_AFs.MGRBphase2final_dpge15_98pct.GnomAD201_dpge15_98pct.GiaB332HC.ccdsexonsplus2bp.csv
