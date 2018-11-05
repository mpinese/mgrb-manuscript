#!/bin/bash
set -euo pipefail

export INDIR=".."
export OUTDIR=".."

mkdir -p tmp/
find "${INDIR}" -maxdepth 1 -name '*.afs.bgz.done' | sed 's/\.afs\.bgz\.done//; s/.*\///g' | sort > tmp/input.txt
find "${OUTDIR}" -maxdepth 1 -name '*.afs.aneu.tsv.done' | sed 's/\.afs\.aneu\.tsv\.done//; s/.*\///g' | sort > tmp/output.txt

function find_aneu
{
  sampleid="${1}"

  echo "${sampleid}"
  ./find-aneuploidy.R \
    --lvt=0.1 \
    --uvt=0.9 \
    --lambdav=5.0 \
    --lambdad=25.0 \
    --madk=10.0 \
    --minrun=5 \
    --dpcalib="../03_average_dp.tsv" \
    --diag="${OUTDIR}/${sampleid}.afs.aneu.diag.pdf" \
    "${INDIR}/${sampleid}.afs.bgz" "${sampleid}" "${OUTDIR}/${sampleid}.afs.aneu.tsv" && \
  touch "${OUTDIR}/${sampleid}.afs.aneu.tsv.done"
}

export -f find_aneu

comm -23 tmp/input.txt tmp/output.txt | parallel find_aneu {}
