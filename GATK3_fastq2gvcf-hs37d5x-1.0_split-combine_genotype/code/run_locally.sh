#!/bin/bash
set -euo pipefail
#set -x

REFERENCE="/nvme/marpin/resources/hs37d5x/reference_genome/hs37d5x.fa"
#REFERENCE="/g/data3/wq2/resources/hs37d5x/reference_genome/hs37d5x.fa"
SHARD_FILE="/nvme/marpin/resources/hs37d5x/reference_genome/hs37d5x.shards"
#SHARD_FILE="/g/data3/wq2/resources/hs37d5x/reference_genome/hs37d5x.shards"
INPUT_DIR="../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine"
OUTPUT_STEM="../local/MGRB.phase2"
TARGET_CHROMS=($(seq 1 22) X Y MT)
TEMP="./tmp"

mkdir -p "../local/"
mkdir -p "${TEMP}"


tier12_all_job_file=$(mktemp "${TEMP}/jobs.tier12.all.XXXXXX")
tier12_match_job_file=$(mktemp "${TEMP}/jobs.tier12.match.XXXXXX")


# Output files:
# MGRB.phase2.tier12. all.     shard${shard_id}.chrom${shard_chrom}.pos${shard_start}-${shard_end}.vcf.gz
# MGRB.phase2.tier123.all.     shard${shard_id}.chrom${shard_chrom}.pos${shard_start}-${shard_end}.vcf.gz
# MGRB.phase2.tier12. sexmatch.shard${shard_id}.chrom${shard_chrom}.pos${shard_start}-${shard_end}.vcf.gz
# MGRB.phase2.tier123.sexmatch.shard${shard_id}.chrom${shard_chrom}.pos${shard_start}-${shard_end}.vcf.gz


function queue_job
{
  shard_id="$1"
  shard_chrom="$2"
  shard_start="$3"
  shard_end="$4"
  input_list_tmpfile="$5"
  outfile="$6"
  jobfile="$7"

  queuefile="${outfile}.queued"
  lockfile="${outfile}.lock"
  donefile="${outfile}.done"

  if [ -e "${queuefile}" ]; then echo "${outfile} already queued"
  elif [ -e "${lockfile}" ]; then echo "${outfile} already running"
  elif [ -e "${donefile}" ]; then echo "${outfile} already done"
  else
    echo -e "${REFERENCE}"\\t"${outfile}"\\t"${shard_chrom}:${shard_start}-${shard_end}"\\t"-V ${input_list_tmpfile}"\\t"gtGVCF.mgrb.qc.${shard_id}" >> "${jobfile}"
#    qsub -z -q normal -vREFERENCE="\"${REFERENCE}\"",OUT_FILE="\"${outfile}\"",SHARD="\"${shard_chrom}:${shard_start}-${shard_end}\"",SOURCE_DESCRIPTOR="\"-V ${input_list_tmpfile}\"" -N "gtGVCF.mgrb.qc.${shard_id}" genotype.pbs
    touch "${queuefile}"
    echo "${outfile} queued"
  fi
}


function queue_tier12_all
{
  shard_id="$1"
  shard_chrom="$2"
  shard_start="$3"
  shard_end="$4"

  outfile="${OUTPUT_STEM}.tier12.all.shard${shard_id}.chrom${shard_chrom}.pos${shard_start}-${shard_end}.vcf.gz"

  input_list_tmpfile=$(mktemp "${TEMP}/gvcf_list.tier12.all.shard${shard_id}.XXXXXX.list")
  ls -1 "${INPUT_DIR}/tier1.match/"*".shard${shard_id}."*.g.vcf.gz > "${input_list_tmpfile}"
  ls -1 "${INPUT_DIR}/tier1.mismatch/"*".shard${shard_id}."*.g.vcf.gz >> "${input_list_tmpfile}"
  ls -1 "${INPUT_DIR}/tier2.match/"*".shard${shard_id}."*.g.vcf.gz >> "${input_list_tmpfile}"
  ls -1 "${INPUT_DIR}/tier2.mismatch/"*".shard${shard_id}."*.g.vcf.gz >> "${input_list_tmpfile}"

  queue_job "${shard_id}" "${shard_chrom}" "${shard_start}" "${shard_end}" "${input_list_tmpfile}" "${outfile}" "${tier12_all_job_file}"
}


function queue_tier12_match
{
  shard_id="$1"
  shard_chrom="$2"
  shard_start="$3"
  shard_end="$4"

  outfile="${OUTPUT_STEM}.tier12.match.shard${shard_id}.chrom${shard_chrom}.pos${shard_start}-${shard_end}.vcf.gz"

  input_list_tmpfile=$(mktemp "${TEMP}/gvcf_list.tier12.match.shard${shard_id}.XXXXXX.list")
  ls -1 "${INPUT_DIR}/tier1.match/"*".shard${shard_id}."*.g.vcf.gz > "${input_list_tmpfile}"
  ls -1 "${INPUT_DIR}/tier2.match/"*".shard${shard_id}."*.g.vcf.gz >> "${input_list_tmpfile}"

  queue_job "${shard_id}" "${shard_chrom}" "${shard_start}" "${shard_end}" "${input_list_tmpfile}" "${outfile}" "${tier12_match_job_file}"
}


# Identify jobs and add to queues
while IFS=$'\t' read -r -a shard_array; do
  shard_chrom="${shard_array[0]}"
  shard_start="${shard_array[5]}"
  shard_end="${shard_array[2]}"
  shard_id="${shard_array[3]}"

  in_target_shards=0
  for target_shard in "$@"; do
    if [ "${shard_id}" == "${target_shard}" ]; then
      in_target_shards=1
      break
    fi
  done

  if [ ${in_target_shards} -eq 1 ]; then
    in_target_chroms=0
    for target_chrom in ${TARGET_CHROMS[@]}; do
      if [ "${shard_chrom}" == "${target_chrom}" ]; then
        in_target_chroms=1
        break
      fi
    done

    if [ ${in_target_chroms} -eq 1 ]; then
#      queue_tier12_all "${shard_id}" "${shard_chrom}" "${shard_start}" "${shard_end}"
      queue_tier12_match "${shard_id}" "${shard_chrom}" "${shard_start}" "${shard_end}"
    fi
  fi
done < "${SHARD_FILE}"


# Execute jobs in queues
if [ $(wc -l "${tier12_match_job_file}" | cut -f1 -d' ') -gt 0 ]; then
  parallel --colsep '\t' ./local_helper.sh \"{1}\" \"{2}\" \"{3}\" "{4}" \"${TEMP}\" :::: "${tier12_match_job_file}"
fi

#if [ $(wc -l "${tier12_all_job_file}" | cut -f1 -d' ') -gt 0 ]; then
#  parallel --colsep '\t' REFERENCE="{1}" OUT_FILE="{2}" SHARD="{3}" SOURCE_DESCRIPTOR="{4}" PBS_JOBFS="${TEMP}" bash genotype.pbs :::: "${tier12_all_job_file}"
#fi

