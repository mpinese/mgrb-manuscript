#!/bin/bash
set -euo pipefail

mkdir -p temp

function frqx2db {
  tmpfile=$(mktemp -p temp/)

  awk '
BEGIN {
  FS="\t"
  OFS=","
}

(NR>1) {
  split($2, varparts, ":")
  if ($3 == varparts[4])
    print varparts[1], varparts[2], varparts[3], varparts[4], $7, $6, $5, $10;
  else
    print varparts[1], varparts[2], varparts[3], varparts[4], $5, $6, $7, $10;
}' \
    < "$1" \
    > "${tmpfile}"

sqlite3 "$2" <<EOF
DROP TABLE IF EXISTS allelefreqs;
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
.mode csv
.import "${tmpfile}" allelefreqs
EOF
}


export -f frqx2db

rm -f temp/02_db_jobs.txt
echo -e '../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/plink/08_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.GiaB_HCR.split.minrep.frqx\t../MGRB.phase2final.GiaB_HCR.split.minrep.db' >> temp/02_db_jobs.txt
echo -e '../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/plink/08b_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.GiaB_HCR.split.minrep.somaticfiltered.frqx\t../MGRB.phase2final.GiaB_HCR.split.minrep.somaticfiltered.db' >> temp/02_db_jobs.txt
echo -e '../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/plink/08_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.m.GiaB_HCR.split.minrep.frqx\t../MGRB.phase2final.m.GiaB_HCR.split.minrep.db' >> temp/02_db_jobs.txt
echo -e '../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/plink/08_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.f.GiaB_HCR.split.minrep.frqx\t../MGRB.phase2final.f.GiaB_HCR.split.minrep.db' >> temp/02_db_jobs.txt
echo -e '../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/plink/07_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.45andup_followup_qcpass_anycancer_mf.GiaB_HCR.split.minrep.frqx\t../45andup.anycancer.mf.GiaB_HCR.split.minrep.db' >> temp/02_db_jobs.txt
echo -e '../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/plink/07_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.45andup_followup_qcpass_breastcancer_f.GiaB_HCR.split.minrep.frqx\t../45andup.breastcancer.f.GiaB_HCR.split.minrep.db' >> temp/02_db_jobs.txt
echo -e '../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/plink/07_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.45andup_followup_qcpass_breastcancer_mf.GiaB_HCR.split.minrep.frqx\t../45andup.breastcancer.mf.GiaB_HCR.split.minrep.db' >> temp/02_db_jobs.txt
echo -e '../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/plink/07_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.45andup_followup_qcpass_colorectalcancer_mf.GiaB_HCR.split.minrep.frqx\t../45andup.colorectalcancer.mf.GiaB_HCR.split.minrep.db' >> temp/02_db_jobs.txt
echo -e '../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/plink/07_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.45andup_followup_qcpass_melanomacancer_mf.GiaB_HCR.split.minrep.frqx\t../45andup.melanomacancer.mf.GiaB_HCR.split.minrep.db' >> temp/02_db_jobs.txt
echo -e '../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/plink/07_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.45andup_followup_qcpass_nocancer_f.GiaB_HCR.split.minrep.frqx\t../45andup.nocancer.f.GiaB_HCR.split.minrep.db' >> temp/02_db_jobs.txt
echo -e '../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/plink/07_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.45andup_followup_qcpass_nocancer_mf.GiaB_HCR.split.minrep.frqx\t../45andup.nocancer.mf.GiaB_HCR.split.minrep.db' >> temp/02_db_jobs.txt
echo -e '../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/plink/07_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.45andup_followup_qcpass_nocancer_m.GiaB_HCR.split.minrep.frqx\t../45andup.nocancer.m.GiaB_HCR.split.minrep.db' >> temp/02_db_jobs.txt
echo -e '../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/plink/07_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.45andup_followup_qcpass_nonmelskincancer_mf.GiaB_HCR.split.minrep.frqx\t../45andup.nonmelskincancer.mf.GiaB_HCR.split.minrep.db' >> temp/02_db_jobs.txt
echo -e '../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/plink/08_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.ch10.GiaB_HCR.split.minrep.frqx\t../MGRB.phase2final.ch10.GiaB_HCR.split.minrep.db' >> temp/02_db_jobs.txt
echo -e '../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/plink/08_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.noch10.GiaB_HCR.split.minrep.frqx\t../MGRB.phase2final.noch10.GiaB_HCR.split.minrep.db' >> temp/02_db_jobs.txt
echo -e '../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/plink/08b_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.ch10.GiaB_HCR.split.minrep.somaticfiltered.frqx\t../MGRB.phase2final.ch10.GiaB_HCR.split.minrep.somaticfiltered.db' >> temp/02_db_jobs.txt
echo -e '../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/plink/08b_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.noch10.GiaB_HCR.split.minrep.somaticfiltered.frqx\t../MGRB.phase2final.noch10.GiaB_HCR.split.minrep.somaticfiltered.db' >> temp/02_db_jobs.txt
echo -e '../../../databases/ASRB/plink/01_ASRB.WGStier12.GiaB_HCR.split.minrep.frqx\t../ASRB.WGStier12.GiaB_HCR.split.minrep.db' >> temp/02_db_jobs.txt
echo -e '../../../databases/ASRB/plink/01_ASRB.WGStier12.GiaB_HCR.split.minrep.m.frqx\t../ASRB.WGStier12.GiaB_HCR.split.minrep.m.db' >> temp/02_db_jobs.txt
echo -e '../../../databases/ASRB/plink/01_ASRB.WGStier12.GiaB_HCR.split.minrep.f.frqx\t../ASRB.WGStier12.GiaB_HCR.split.minrep.f.db' >> temp/02_db_jobs.txt
echo -e '../../../databases/ASRB/plink/01b_ASRB.WGStier12.GiaB_HCR.split.minrep.somaticfiltered.frqx\t../ASRB.WGStier12.GiaB_HCR.split.minrep.somaticfiltered.db' >> temp/02_db_jobs.txt
echo -e '../../../databases/ASRB/plink/01b_ASRB.WGStier12.GiaB_HCR.split.minrep.somaticfiltered.m.frqx\t../ASRB.WGStier12.GiaB_HCR.split.minrep.somaticfiltered.m.db' >> temp/02_db_jobs.txt
echo -e '../../../databases/ASRB/plink/01b_ASRB.WGStier12.GiaB_HCR.split.minrep.somaticfiltered.f.frqx\t../ASRB.WGStier12.GiaB_HCR.split.minrep.somaticfiltered.f.db' >> temp/02_db_jobs.txt


parallel --colsep '\t' frqx2db {1} {2} :::: temp/02_db_jobs.txt
