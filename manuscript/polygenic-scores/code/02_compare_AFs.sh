#!/bin/bash
set -euo pipefail

python3 02h_generate_MGRB_AF_db.py

bash ./02h_generate_MGRB_GnomAD_AF_db_csv.sh

Rscript 02h_clean_ebi_gwas_data.R

mkdir -p log
Rscript 02h_test_AF_differences.R 2>&1 | tee log/02h_test_AF_differences.log
