#!/bin/bash
set -euo pipefail

export INDIR="../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.commonhets"
export OUTDIR=".."

mkdir -p tmp/
find "${INDIR}" -maxdepth 1 -name '*.afs.bgz.done' | sed 's/\.afs\.bgz\.done//; s/.*\///g' | sort > tmp/input.txt
find "${OUTDIR}" -maxdepth 1 -name '*.afs.aneu2.rds.done' | sed 's/\.afs\.aneu2\.rds\.done//; s/.*\///g' | sort > tmp/output.txt

function find_aneu
{
  sampleid="${1}"

  echo "${sampleid}"
  ./find-aneuploidy2.R \
    "../03_average_dp.tsv" \
    "../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.commonhets.gc.tsv" \
    "${OUTDIR}/${sampleid}.afs.aneu2.diag.pdf" \
    "${INDIR}/${sampleid}.afs.bgz" \
    "${OUTDIR}/${sampleid}.afs.aneu2.rds" && \
  touch "${OUTDIR}/${sampleid}.afs.aneu2.rds.done"
}

export -f find_aneu

comm -23 tmp/input.txt tmp/output.txt | parallel find_aneu {}

