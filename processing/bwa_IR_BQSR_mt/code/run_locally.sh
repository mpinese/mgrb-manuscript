#!/bin/bash
set -euo pipefail

INDIR="../../GATK3_fastq2gvcf-hs37d5x-1.0/bams"
OUTDIR=".."


function execute {
  chrom=$2
  infile="${1%.md5.OK}"
  filebase=$(basename "${infile}")
  sampleid="${filebase%%_*}"

  output_file="${OUTDIR}/${filebase%.bam}.${chrom}.bam"

  lock_file="${output_file}.lock"
  done_file="${output_file}.done"
  log_file="${output_file}.log"

  if [ -e "${done_file}" ]; then
    echo "${sampleid} ${chrom} already done."
    return
  elif [ -e "${lock_file}" ]; then
    echo "${sampleid} ${chrom} locked."
    return
  fi

  touch "${lock_file}"

  ./bin/samtools view -b -o "${output_file}" "${infile}" ${chrom} 2>&1 | tee -a "${log_file}"
  ./bin/samtools index "${output_file}" 2>&1 | tee -a "${log_file}"

  touch "${done_file}"
  rm -f "${lock_file}"

  echo "${sampleid} ${chrom} done"
}

export -f execute
export INDIR
export OUTDIR


#parallel execute {} MT ::: "${INDIR}"/*.bam.md5.OK
parallel execute {} NC_001422.1 ::: "${INDIR}"/*.bam.md5.OK

