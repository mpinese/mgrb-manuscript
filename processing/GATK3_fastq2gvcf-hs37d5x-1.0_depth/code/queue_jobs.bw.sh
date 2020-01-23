#!/bin/bash
set -euo pipefail

INDIR="../../GATK3_fastq2gvcf-hs37d5x-1.0/bams"
OUTDIR=".."


for md5file in "${INDIR}"/*.bam.md5.OK; do
    file="${md5file%.md5.OK}"
    sampleid=$(basename "${file}")
    sampleid="${sampleid%%_*}"

    output_file_starch="${OUTDIR}/${sampleid}.depth.starch"

    queue_file="${OUTDIR}/${sampleid}.depth.starch.queued"
    lock_file="${OUTDIR}/${sampleid}.depth.starch.lock"
    done_file="${OUTDIR}/${sampleid}.depth.starch.done"
    term_file="${OUTDIR}/${sampleid}.depth.starch.term"
    log_file="${OUTDIR}/${sampleid}.depth.starch.log"

    if [ -e "${queue_file}" ]; then
        echo "${sampleid} already queued"
    elif [ -e "${lock_file}" ]; then
        echo "${sampleid} already running"
    elif [ -e "${done_file}" ]; then
        echo "${sampleid} already done"
    elif [ -e "${term_file}" ]; then
        echo "${sampleid} was terminated"
    else
        qsub -q expressbw -z -v INFILE="${file}",OUTFILE="${output_file_starch}",SAMPLEID="${sampleid}" -N dp${sampleid} samtools-depth.bw.pbs
        touch "${queue_file}"
        echo "${sampleid} queued"
    fi
done
