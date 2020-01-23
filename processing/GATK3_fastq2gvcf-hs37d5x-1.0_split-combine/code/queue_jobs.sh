#!/bin/bash
set -euo pipefail

mkdir -p ../tier1.match/
mkdir -p ../tier2.match/
mkdir -p ../tier3.match/
mkdir -p ../tier1.mismatch/
mkdir -p ../tier2.mismatch/
mkdir -p ../tier3.mismatch/
mkdir -p ../tier4.all/

bash queue_jobs.byshard.bygroup.bylist.sh ../tier1.match/MGRB.phase2.tier1.match ../../GATK3_fastq2gvcf-hs37d5x-1.0_qc_split-combine_genotype_sampleqc/tier1.sexmatch.sample_list 100
bash queue_jobs.byshard.bygroup.bylist.sh ../tier2.match/MGRB.phase2.tier2.match ../../GATK3_fastq2gvcf-hs37d5x-1.0_qc_split-combine_genotype_sampleqc/tier2.sexmatch.sample_list 100
bash queue_jobs.byshard.bygroup.bylist.sh ../tier3.match/MGRB.phase2.tier3.match ../../GATK3_fastq2gvcf-hs37d5x-1.0_qc_split-combine_genotype_sampleqc/tier3.sexmatch.sample_list 100
bash queue_jobs.byshard.bygroup.bylist.sh ../tier1.mismatch/MGRB.phase2.tier1.mismatch ../../GATK3_fastq2gvcf-hs37d5x-1.0_qc_split-combine_genotype_sampleqc/tier1.sexmismatch.sample_list 100
bash queue_jobs.byshard.bygroup.bylist.sh ../tier2.mismatch/MGRB.phase2.tier2.mismatch ../../GATK3_fastq2gvcf-hs37d5x-1.0_qc_split-combine_genotype_sampleqc/tier2.sexmismatch.sample_list 100
bash queue_jobs.byshard.bygroup.bylist.sh ../tier3.mismatch/MGRB.phase2.tier3.mismatch ../../GATK3_fastq2gvcf-hs37d5x-1.0_qc_split-combine_genotype_sampleqc/tier3.sexmismatch.sample_list 100
bash queue_jobs.byshard.bygroup.bylist.sh ../tier4.all/MGRB.phase2.tier4.all ../../GATK3_fastq2gvcf-hs37d5x-1.0_qc_split-combine_genotype_sampleqc/tier4.all.sample_list 100

