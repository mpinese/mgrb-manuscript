#!/bin/bash
set -euo pipefail

REFERENCE="/g/data3/wq2/resources/reference_genomes/hs37d5x/hs37d5x.fa"
INDIR="../../GATK3_fastq2gvcf-hs37d5x-1.0/bams"
OUTDIR=".."


for md5file in "${INDIR}"/*.bam.md5.OK; do
    file="${md5file%.md5.OK}"
    sampleid=$(basename "${file}")
    sampleid="${sampleid%%_*}"

    output_file_vcfgz="${OUTDIR}/${sampleid}.gridss.vcf.gz"
    output_file_bam="${OUTDIR}/${sampleid}.gridss.bam"

    queue_file="${OUTDIR}/${sampleid}.gridss.vcf.gz.queued"
    lock_file="${OUTDIR}/${sampleid}.gridss.vcf.gz.lock"
    done_file="${OUTDIR}/${sampleid}.gridss.vcf.gz.done"
    term_file="${OUTDIR}/${sampleid}.gridss.vcf.gz.term"
    log_file="${OUTDIR}/${sampleid}.gridss.vcf.gz.log"

    if [ -e "${queue_file}" ]; then
        echo "${sampleid} already queued"
    elif [ -e "${lock_file}" ]; then
        echo "${sampleid} already running"
    elif [ -e "${done_file}" ]; then
        echo "${sampleid} already done"
    elif [ -e "${term_file}" ]; then
        echo "${sampleid} was terminated"
    else
        qsub -z -v INFILE="${file}",OUTFILE_VCFGZ="${output_file_vcfgz}",OUTFILE_BAM="${output_file_bam}",SAMPLEID="${sampleid}",REFERENCE="${REFERENCE}" -N gr${sampleid} gridss.pbs
        touch "${queue_file}"
        echo "${sampleid} queued"
    fi
done
