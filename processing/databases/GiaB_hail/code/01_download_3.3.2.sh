#!/bin/bash
set -euo pipefail

mkdir -p source_data
cd source_data
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/README_NISTv3.3.2.txt
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz.tbi
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed

gzip -dc HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz | \
    awk 'BEGIN{FS="\t";OFS="\t"} /^#/; /^[^#]/ { $8="."; n_format_parts=split($9, format_parts, ":"); n_data_parts=split($10, data_parts, ":"); for (i = 1; i <= n_format_parts; i++) { if (format_parts[i] == "GT") { break } }; if (i == n_format_parts + 1) { next }; $9 = format_parts[i]; gt = data_parts[i]; gsub(/\|/, "/", gt); $10=gt; n_gt_parts = split(gt, gt_parts, "/"); for (i = 1; i <= n_gt_parts; i++) { gt_parts[i] = int(gt_parts[i]) }; asort(gt_parts); $10 = gt_parts[1]; for (i = 2; i <= n_gt_parts; i++) { $10 = $10 "/" gt_parts[i] }; print }' | \
    grep -v '^##INFO' | \
    bgzip > HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.minimal.vcf.bgz
tabix HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.minimal.vcf.bgz

cd ..
