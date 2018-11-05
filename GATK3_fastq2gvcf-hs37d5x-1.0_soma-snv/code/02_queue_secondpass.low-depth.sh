#!/bin/bash
set -euo pipefail

# Valid settings found:
# QUEUE=normalbw
# CORES=28
# RAMPERCORE=6
# WALLTIME="1:00:00"
# MAXQUEUED=100

# QUEUE=normalbw
# CORES=14
# RAMPERCORE=8
# WALLTIME="1:45:00"
# MAXQUEUED=100

# QUEUE=normal
# CORES=16
# RAMPERCORE=8
# WALLTIME=???
# MAXQUEUED=100

QUEUE=normalbw
CORES=14
RAMPERCORE=8
WALLTIME="1:45:00"
MAXQUEUED=100
INDIR="../../GATK3_fastq2gvcf-hs37d5x-1.0/bams"
OUTDIR=".."
SHARDS="../../../../../resources/reference_genomes/hs37d5x/hs37d5x.shards"
REFERENCE="../../../../../resources/reference_genomes/hs37d5x/hs37d5x.fa"
BLACKLIST="loh-pass2-commonsnps-blacklist.bed"

num_queued=$((MAXQUEUED+1))

find "${INDIR}" -maxdepth 1 -type f -name '*.bam' | while read input_bam; do
  if [ ${num_queued} -gt ${MAXQUEUED} ]; then
    set +e
    num_queued=$(qstat -u `whoami` | grep normalbw | awk '($10 == "Q"){total+=1} END{print total+0}')
    set -e
  fi
  if [ ${num_queued} -gt ${MAXQUEUED} ]; then
    break
  fi

  sample_name=$(basename "${input_bam}")
  sample_name="${sample_name%%_*}"

  mkdir -p "${OUTDIR}/${sample_name}/pass2/shards/"

  outfile_stem="${OUTDIR}/${sample_name}/pass2/shards/${sample_name}.soma-snv.pass2"

  queuefile="${outfile_stem}.queued"
  lockfile="${outfile_stem}.lock"
  termfile="${outfile_stem}.term"
  donefile="${outfile_stem}.done"

  if [ -e "${queuefile}" ]; then
    echo -e "\e[94m${sample_name} already queued\e[39m"
  elif [ -e "${lockfile}" ]; then
    echo -e "\e[33m${sample_name} running\e[39m"
  elif [ -e "${termfile}" ]; then
    echo -e "\e[91m${sample_name} terminated\e[39m"
  elif [ -e "${donefile}" ]; then
    echo -e "\e[92m${sample_name} already done\e[39m"
  else
    qsub -z -q ${QUEUE} -l walltime="${WALLTIME}" -N "s2.${sample_name}" -l ncpus=${CORES} -l mem=$((RAMPERCORE*CORES))G -v CORES="${CORES}",INFILE="${input_bam}",OUTFILE_STEM="${outfile_stem}",SHARDS="${SHARDS}",BLACKLIST="${BLACKLIST}",REFERENCE="${REFERENCE}" soma-snv.low-depth.pbs && \
    echo -e "\e[39m${sample_name} queued\e[39m" && \
    touch "${queuefile}" && \
    num_queued=$((num_queued+1)) && \
    sleep 0.5
  fi
done

