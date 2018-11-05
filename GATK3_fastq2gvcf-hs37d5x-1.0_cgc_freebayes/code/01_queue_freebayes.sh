#!/bin/bash
set -euo pipefail

#QUEUE=normalbw
QUEUE=expressbw
MAXQUEUED=100

INDIR="../../GATK3_fastq2gvcf-hs37d5x-1.0/bams"
OUTDIR=".."

REFERENCE="../../../../../resources/hs37d5x/reference_genome/hs37d5x.fa"
TARGETS="cosmic_cgc_20171211.refgene.pm10kb.bed"

num_queued=$((MAXQUEUED+1))

find "${INDIR}" -maxdepth 1 -name '*.bam' | while read inbam; do
  if [ ${num_queued} -gt ${MAXQUEUED} ]; then
    set +e
    num_queued=$(qstat -u `whoami` | grep "${QUEUE:0:8}" | awk '($10 == "Q"){total+=1} END{print total+0}')
    set -e
  fi
  if [ ${num_queued} -gt ${MAXQUEUED} ]; then
    echo -e "\e[31mTerminating: queue full\e[39m"
    break
  fi

  fname=$(basename "${inbam}")
  sname="${fname%%_*}"
  outvcfgz="${OUTDIR}/${fname%.bam}.cgc.freebayes.vcf.gz"

  queue_file="${outvcfgz}.queued"
  done_file="${outvcfgz}.done"
  lock_file="${outvcfgz}.lock"

  if [ -e "${queue_file}" ]; then
    echo "${outvcfgz} already queued"
  elif [ -e "${done_file}" ]; then
    echo "${outvcfgz} already done"
  elif [ -e "${lock_file}" ]; then
    echo "${outvcfgz} already running"
  else
    qsub -z -N "fb${sname}" -q "${QUEUE}" -v REFERENCE="${REFERENCE}",INBAM="${inbam}",INBED="${TARGETS}",OUTVCFGZ="${outvcfgz}",QUEUEFILE="${queue_file}",LOCKFILE="${lock_file}",DONEFILE="${done_file}" freebayes.pbs && \
    echo "${outvcfgz} queued" && \
    touch "${queue_file}" && \
    num_queued=$((num_queued+1)) && \
    sleep 1
  fi
done


