#!/bin/bash
set -euo pipefail

DEPTH_STATS_DIR="../../GATK3_fastq2gvcf-hs37d5x-1.0_depth_stats"
SAMPLEID_DIR="../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail"

# Requires:
#  Good sample IDs in ../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/tier[12].sample_list
#  Depth stats files in ../../GATK3_fastq2gvcf-hs37d5x-1.0_depth_stats/

required_paths=("${SAMPLEID_DIR}/tier1.sample_list" "${SAMPLEID_DIR}/tier2.sample_list" "${DEPTH_STATS_DIR}")

for required_path in ${required_paths[@]}; do
  if [ ! -e "${required_path}" ]; then
    echo -e "\e[31mCannot find required path ${required_path}\e[0m"
    exit 1
  fi
done

# Concatenate all available <sampleID>.dpgood.bed.gz files for good samples
mkdir -p tmp/
rm -f tmp/dpgood_paths.txt
good_depth_sample_count=0
cat "${SAMPLEID_DIR}/tier1.sample_list" "${SAMPLEID_DIR}/tier2.sample_list" | while read good_sample_id; do
  good_sample_dpgood_path="${DEPTH_STATS_DIR}/${good_sample_id}.dpgood.bed.gz"
  if [ -e "${good_sample_dpgood_path}" ]; then
    echo "${good_sample_dpgood_path}" >> tmp/dpgood_paths.txt
    good_depth_sample_count=$((good_depth_sample_count+1))
  else
    echo -e "\e[93mWarning: missing dpgood file for sample ${good_sample_id}\e[0m" > /dev/stderr
  fi
done

echo -e "\e[36mExtracting dpgood files...\e[0m"
mkdir -p tmp/0
rm -f tmp/0/*
parallel "gzip -dc {} > tmp/0/{/.}" :::: tmp/dpgood_paths.txt

echo -e "\e[36mCounting overlaps...\e[0m"
find tmp/0 -name '*.bed' -print0 | LC_ALL=C sort -T tmp/ -k1,1 -k2,2n -m --files0-from='-' --batch-size 60 | bedtools genomecov -i - -g <(gzip -dc ../hs37d5x_data/genome.bed.gz | cut -f 1,3) -bga | bgzip > ../gooddepth_samplecounts.bed.gz
wc -l tmp/dpgood_paths.txt | cut -d' ' -f1 > ../gooddepth_totalsamples.txt

echo -e "\e[36mDefining bad regions...\e[0m"
# frac(good_depth) < 0.95 --> morphological close (size 65) --> morphological open (size 5)
#
# Justification: we wish to exclude regions of possible false
# positive variation due to cryptic mismapping.
# Such regions appear as:
#
# Mismapping fragment:
#
#   @@@@++++++@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@+++++@@
#
#
# Read depth:
#
#       High dp                               High dp
#         |                                     |
#         V                                     V
#
#       ######                                #####
#       ######                                #####
#       ######                                #####
# #####################################################
#
# Reference:
# -----------------------------------------------------
#
# Two regions of high depth are present, corresponding
# to regions (+) on a mismappping fragment (@) that
# coincidentally map perfectly to the reference (-).
# The concern is that a read that captures both ends
# of the mismapping fragment will be given a high enough
# MQ, and low enough soft clipped fraction, to be considered
# by the caller, perhaps as an indel.  Such reads will
# span both high dp regions (else it will be a terminal
# soft clipped read, and likely ignored by the caller).
# For a read length of 150 bp, and a high dp matching
# region minimum size of 10 bp, this corresponds to a
# maximum spanning space size of (150-10-10) = 130 bp.
# For a morphological closing operation which closes
# from both sides, the filter size is then 130/2 = 65 bp.
gzip -dc ../gooddepth_samplecounts.bed.gz | \
  awk "BEGIN{FS=\"\t\";OFS=\"\t\"} (\$4 / $(cat ../gooddepth_totalsamples.txt) < 0.95);" | \
  awk 'BEGIN{FS="\t";OFS="\t"} {$2 -= 65; $3 += 65; if ($2 < 0) {$2 = 0}; print $1,$2,$3}' | \
  bedtools merge | \
  awk 'BEGIN{FS="\t";OFS="\t"} {print $1,$2+65,$3-65}' | \
  awk 'BEGIN{FS="\t";OFS="\t"} ($3-$2 > 10);' | \
  bgzip > ../hs37d5x_data/bad.depth.bed.gz

