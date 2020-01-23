#!/bin/bash
set -euo pipefail

REFERENCE="/g/data3/wq2/resources/hs37d5x/reference_genome/hs37d5x.fa"
SHARD_FILE="/g/data3/wq2/resources/hs37d5x/qc_loci/qc_loci.InfiniumQCArray-24v1-0_A3.shards"
LOCUS_FILE="/g/data3/wq2/resources/hs37d5x/qc_loci/qc_loci.InfiniumQCArray-24v1-0_A3.bed"
INPUT_DIR="/g/data3/wq2/results/phase2/hs37d5x/GATK3_fastq2gvcf-hs37d5x-1.0_qc_split-combine"
OUTPUT_STEM="../MGRB.phase2qc"
TARGET_CHROMS=($(seq 1 22) X Y MT)
TEMP="./tmp"

mkdir -p "${TEMP}"

while IFS=$'\t' read -r -a shard_array; do
  shard_chrom="${shard_array[0]}"
  shard_start="${shard_array[5]}"
  shard_end="${shard_array[2]}"
  shard_id="${shard_array[3]}"

  shard_tmpfile=$(mktemp "${TEMP}/shard${shard_id}.XXXXXX.bed")
  echo -e "${shard_chrom}\t$((shard_start-1))\t${shard_end}\n" > "${shard_tmpfile}"
  shard_locus_tmpfile=$(mktemp "${TEMP}/locus_list.shard${shard_id}.XXXXXX.bed")
  bedtools intersect -wa -a "${LOCUS_FILE}" -b "${shard_tmpfile}" > "${shard_locus_tmpfile}"

  in_target_chroms=0
  for target_chrom in ${TARGET_CHROMS[@]}; do
    if [ "${shard_chrom}" == "${target_chrom}" ]; then
      in_target_chroms=1
      break
    fi
  done

  if [ ${in_target_chroms} -eq 1 ]; then
    outfile="${OUTPUT_STEM}.shard${shard_id}.chrom${shard_chrom}.pos${shard_start}-${shard_end}.vcf.gz"
    queuefile="${outfile}.queued"
    lockfile="${outfile}.lock"
    donefile="${outfile}.done"

    input_list_tmpfile=$(mktemp "${TEMP}/gvcf_list.shard${shard_id}.XXXXXX.list")
    ls -1 "${INPUT_DIR}/"*".shard${shard_id}."*.g.vcf.gz > "${input_list_tmpfile}"

    if [ -e "${queuefile}" ]; then echo "${outfile} already queued"
    elif [ -e "${lockfile}" ]; then echo "${outfile} already running"
    elif [ -e "${donefile}" ]; then echo "${outfile} already done"
    else
      qsub -z -q express -vREFERENCE="\"${REFERENCE}\"",OUT_FILE="\"${outfile}\"",SHARD="\"${shard_locus_tmpfile}\"",SOURCE_DESCRIPTOR="\"-V ${input_list_tmpfile}\"" -N "gtGVCF.mgrb.qc.${shard_id}" genotype.pbs
      touch "${queuefile}"
      echo "${outfile} queued"
    fi
  fi
done < "${SHARD_FILE}"
