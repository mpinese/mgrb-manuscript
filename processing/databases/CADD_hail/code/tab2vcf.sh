#!/bin/bash
set -euo pipefail

infile="$1"
outfile="$2"
regionfile="$3"
interval="$4"

mkdir -p tmp

tempfile1=$(mktemp tmp/XXXXXXXXX.vcf.gz)
tempfile2=$(mktemp tmp/XXXXXXXXX.vcf.gz)

tabix "${regionfile}" "${interval}" | while IFS=$'\t' read -r chrom start0 stop1; do
  tabix "${infile}" "${chrom}:$((start0+1))-${stop1}"
done | \
  awk '\
    BEGIN {
      FS="\t"
      OFS="\t"
      print "##fileformat=VCFv4.2"
      print "##INFO=<ID=CADDRaw,Number=1,Type=Float,Description=\"CADD raw score\">"
      print "##INFO=<ID=CADDPhred,Number=1,Type=Float,Description=\"CADD Phred-scaled score\">"
      print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    }

    /^[^#]/ {
      print $1,$2,".",$3,$4,".",".","CADDRaw="$5";CADDPhred="$6
    }' | \
  bgzip > "${tempfile1}"

tabix "${tempfile1}"

./bin/vt normalize -r ../../../resources/hs37d5x/reference_genome/hs37d5x.fa -o "${tempfile2}" "${tempfile1}"

mv "${tempfile2}" "${outfile}"
tabix "${outfile}"


