#!/bin/bash
set -euo pipefail

# Called as
# bash queue_jobs.byshard.bygroup.bylist.sh <outputstem> <idlist> <batchsize>

OUTPUT_STEM="$1"
SAMPLEID_LIST="$2"
BATCH_SIZE="$3"

REFERENCE="/g/data3/wq2/resources/hs37d5x/reference_genome/hs37d5x.fa"
SHARD_FILE="/g/data3/wq2/resources/hs37d5x/reference_genome/hs37d5x.shards"
INPUT_DIR="/g/data3/wq2/results/phase2/hs37d5x/GATK3_fastq2gvcf-hs37d5x-1.0/gvcfs"
TARGET_CHROMS=($(seq 1 22) X Y MT)
TEMP="./tmp"

mkdir -p "${TEMP}"

function queue_batch
{
  sources="$1"
  sid1="$2"
  sid2="$3"

  while IFS=$'\t' read -r -a shard_array; do
    shard_chrom="${shard_array[0]}"
    shard_start="${shard_array[5]}"
    shard_end="${shard_array[2]}"
    shard_id="${shard_array[3]}"

    in_target_chroms=0
    for target_chrom in ${TARGET_CHROMS[@]}; do
      if [ "${shard_chrom}" == "${target_chrom}" ]; then
        in_target_chroms=1
        break
      fi
    done

    if [ ${in_target_chroms} -eq 1 ]; then
      outfile="${OUTPUT_STEM}.group${sid1}-${sid2}.shard${shard_id}.chrom${shard_chrom}.pos${shard_start}-${shard_end}.g.vcf.gz"
      queuefile="${outfile}.queued"
      lockfile="${outfile}.lock"
      donefile="${outfile}.done"
      termfile="${outfile}.term"

      if [ -e "${queuefile}" ]; then echo "  ${outfile} already queued"
      elif [ -e "${lockfile}" ]; then echo "  ${outfile} already running"
      elif [ -e "${donefile}" ]; then echo "  ${outfile} already done"
      elif [ -e "${termfile}" ]; then echo "  ${outfile} terminated"
      else
        qsub -z -vREFERENCE="\"${REFERENCE}\"",OUT_FILE="\"${outfile}\"",SHARD="\"${shard_chrom}:${shard_start}-${shard_end}\"",SOURCE_DESCRIPTOR="\"${sources}\"" -N SCgVCF split_combine.pbs
        touch "${queuefile}"
        echo "  ${outfile} queued"
      fi
    fi
  done < "${SHARD_FILE}"
}

i=0
source_descriptor=""
sampleid=""
while read line; do
  if [ -z "${line// }" ]; then
    continue
  fi
  sampleid="${line}"

  echo "  $((i+1))/$BATCH_SIZE ${sampleid}"
  if [ "$i" -eq 0 ]; then
    first_sid="${sampleid}"
  fi

  if [[ "${source_descriptor}" == "" ]]; then
    input_list_tmpfile=$(mktemp "${TEMP}/vcf_list.group${sampleid}-.XXXXXX.list")
    source_descriptor="-V \"${input_list_tmpfile}\""
  fi

  echo "${INPUT_DIR}/${sampleid}_"*".g.vcf.gz" >> "${input_list_tmpfile}"
  i=$((i+1))

  if [ "$i" -eq "${BATCH_SIZE}" ]; then
    echo "Sending jobs..."
    queue_batch "${source_descriptor}" "${first_sid}" "${sampleid}"
    i=0
    source_descriptor=""
  fi
done < "${SAMPLEID_LIST}"

if [ "$i" -gt 0 ]; then
  queue_batch "${source_descriptor}" "${first_sid}" "${sampleid}"
fi
