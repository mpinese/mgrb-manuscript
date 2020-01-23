#!/bin/bash
set -euo pipefail

PROJECT=wq2
QUEUE=normalbw
PARTASKS=28
NCORES=28
WALLTIME=6:00:00

INDIR="../../GATK3_fastq2gvcf-hs37d5x-1.0/bams/"
OUTDIR=".."
GATK="bin/GenomeAnalysisTK.jar"
REFERENCE="../../../../../resources/hs37d5x/reference_genome/hs37d5x.fa"
LOCI="data/soma-loh.hs37.loci.tsv.xz"
GC="data/soma-loh.hs37.gc.tsv.xz"
AFFINITY="data/soma-loh.hs37.affinity.mgrb.tsv.xz"


mkdir -p tmp

find "${INDIR}" -maxdepth 1 -name '*.bam' > tmp/paths_avail.txt
while read path; do basename "${path}"; done < tmp/paths_avail.txt | sed 's/\..*//' | sort > tmp/samples_avail.txt
find "${OUTDIR}" -maxdepth 1 -name '.done' | while read path; do basename "${path}"; done | sed 's/\..*//' | sort > tmp/samples_done.txt
comm -23 tmp/samples_avail.txt tmp/samples_done.txt > tmp/samples_todo.txt
grep -Ff tmp/samples_todo.txt tmp/paths_avail.txt > tmp/paths_todo.txt
while read path; do
  sampleid=$(basename "${path}")
  sampleid="${sampleid%%.dupmarked.*}"
  sampleid="${sampleid//\./_}"
  if [ -e "${OUTDIR}/${sampleid}.done" ] | [ -e "${OUTDIR}/${sampleid}.lock" ] | [ -e "${OUTDIR}/${sampleid}.queued" ]; then
    continue
  fi
  echo -e "${path}\t${OUTDIR}/${sampleid}"
done < tmp/paths_todo.txt > tmp/jobs_todo.txt

mkdir -p tmp/shard$$
split -l 56 tmp/jobs_todo.txt tmp/shard$$/

for f in tmp/shard$$/*; do
  qsub -z -P wq2 -q "${QUEUE}" -l walltime="${WALLTIME}" -l wd -l other=gdata3 -l mem=128G -l ncpus=${NCORES} -N sl -v GATK="${GATK}",LOCI="${LOCI}",GC="${GC}",AFFINITY="${AFFINITY}",REFERENCE="${REFERENCE}",JOBS="${f}",PARTASKS="${PARTASKS}" soma-loh.pbs && \
  { cut -f 2 "${f}" | while read g; do touch "${g}.queued"; done } && \
  echo "${f} queued" && \
  sleep 1
done

