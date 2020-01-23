#!/bin/bash
set -euo pipefail

REFERENCE="/g/data3/wq2/resources/reference_genomes/hs37d5x/hs37d5x.fa"
INDIR="../../GATK3_fastq2gvcf-hs37d5x-1.0/bams"
OUTDIR=".."


for md5file in "${INDIR}"/*.bam.md5.OK; do
    file="${md5file%.md5.OK}"
    filebase=$(basename "${file}")
    sampleid="${filebase%%_*}"

    output_file="${OUTDIR}/${filebase%.bam}.MT.bam"
    queue_file="${output_file}.queued"
    lock_file="${output_file}.lock"
    done_file="${output_file}.done"
    term_file="${output_file}.term"
    log_file="${output_file}.log"
    if [ -e "${queue_file}" ]; then
        echo "${sampleid} MT already queued"
    elif [ -e "${lock_file}" ]; then
        echo "${sampleid} MT already running"
    elif [ -e "${done_file}" ]; then
        echo "${sampleid} MT already done"
    elif [ -e "${term_file}" ]; then
        echo "${sampleid} MT was terminated"
    else
        qsub -z -v INFILE="${file}",OUTFILE="${output_file}",CHROMS=MT -N mt${sampleid} samtools_view.pbs
        touch "${queue_file}"
        echo "${sampleid} MT queued"
    fi

    output_file="${OUTDIR}/${filebase%.bam}.phiX.bam"
    queue_file="${output_file}.queued"
    lock_file="${output_file}.lock"
    done_file="${output_file}.done"
    term_file="${output_file}.term"
    log_file="${output_file}.log"
    if [ -e "${queue_file}" ]; then
        echo "${sampleid} PhiX already queued"
    elif [ -e "${lock_file}" ]; then
        echo "${sampleid} PhiX already running"
    elif [ -e "${done_file}" ]; then
        echo "${sampleid} PhiX already done"
    elif [ -e "${term_file}" ]; then
        echo "${sampleid} PhiX was terminated"
    else
        qsub -z -v INFILE="${file}",OUTFILE="${output_file}",CHROMS=NC_001422.1 -N px${sampleid} samtools_view.pbs
        touch "${queue_file}"
        echo "${sampleid} PhiX queued"
    fi
done

