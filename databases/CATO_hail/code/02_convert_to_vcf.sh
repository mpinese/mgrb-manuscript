#!/bin/bash
set -euo pipefail

mkdir -p ../converted_data
unstarch ../source_data/dbSNP142.CATO.V1.1.starch | awk '\
  BEGIN {
    FS="\t"
    OFS="\t"
    print "##fileformat=VCFv4.2"
    print "##source=dbSNP142.CATO.V1.1.txt.gz"
    print "##INFO=<ID=score,Number=1,Type=Float,Description=\"CATO raw score\">"
    print "##INFO=<ID=motif,Number=1,Type=String,Description=\"Overlapping motif\">"
    print "##INFO=<ID=strand,Number=1,Type=Integer,Description=\"Strand on overlapping motif\">"
    print "##INFO=<ID=motifpos,Number=1,Type=Integer,Description=\"Position within overlapping motif\">"
    print "##INFO=<ID=celltypes,Number=.,Type=String,Description=\"Cell types for prediction\">"
    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  }

  (NR > 1) {
    chrom=$1
    sub(/^chr/, "", chrom)

    strand=0
    if ($6 == "+")
        strand=+1
    else if ($6 == "-")
        strand=-1

    celltypes=$11
    gsub(/;/, ",", celltypes)

    print chrom,$3,$4,$9,$10,".",".","score="$5";motif="$7";strand="strand";motifpos="$8";celltypes="celltypes
  }' | bgzip > ../converted_data/dbSNP142.CATO.V1.1.vcf.gz

tabix ../converted_data/dbSNP142.CATO.V1.1.vcf.gz

./bin/vt normalize -r ../../../resources/hs37d5x/reference_genome/hs37d5x.fa -o ../converted_data/dbSNP142.CATO.V1.1.norm.vcf.gz ../converted_data/dbSNP142.CATO.V1.1.vcf.gz

mv ../converted_data/dbSNP142.CATO.V1.1.norm.vcf.gz ../converted_data/dbSNP142.CATO.V1.1.norm.vcf.bgz
tabix ../converted_data/dbSNP142.CATO.V1.1.norm.vcf.bgz
