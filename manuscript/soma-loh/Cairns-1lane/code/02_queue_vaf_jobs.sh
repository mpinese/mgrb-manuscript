#!/bin/bash
set -euo pipefail
#set -x

INDIR="../../bwa"
OUTDIR=".."

for input_donefile in "${INDIR}"/*.bam.done; do
  infile="${input_donefile%.done}"
  outfile=$(basename "${infile}")
  outfile="${OUTDIR}/${outfile%.bam}.afs.bgz"

  donefile="${outfile}.done"
  termfile="${outfile}.term"
  lockfile="${outfile}.lock"
  queuefile="${outfile}.queued"

  if [ -e "${donefile}" ]; then
    echo "${outfile} already done."
  elif [ -e "${termfile}" ]; then
    echo "${outfile} terminated."
  elif [ -e "${lockfile}" ]; then
    echo "${outfile} running."
  elif [ -e "${queuefile}" ]; then
    echo "${outfile} already queued."
  else
    qsub -z -v INFILE="${infile}",OUTFILE="${outfile}" vaf_worker.pbs && \
    touch "${queuefile}" && \
    echo "${outfile} queued" && \
    sleep 0.5
  fi
done

