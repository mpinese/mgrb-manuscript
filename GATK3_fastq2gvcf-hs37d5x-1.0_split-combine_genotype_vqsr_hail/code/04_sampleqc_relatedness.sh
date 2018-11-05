#!/bin/bash
set -euo pipefail

./bin/pyhail.sh ./04h_ldprune_for_pcrelate.py

# Run PCRELATE (GENESIS package) on the high-quality SNP PLINK BED generated in 04_merge_with_ref_cohorts.py
Rscript --slave ./04h_pcrelate.R | tee log/04h_pcrelate.R.log
# Generates various data and plots, but primarily the set of samples to drop to remove 2nd-degree
# relatives and closer:
# ../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.1000G_HRC_commonSNVs.ldpruned.pcrelate.samples_to_drop

./bin/pyhail.sh ./04h_sampleqc_drop_related.py
