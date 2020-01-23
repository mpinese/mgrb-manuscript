#!/bin/bash
set -euo pipefail

mkdir -p ../converted_data
gzip -dc ../source_data/Eigen_hg19_coding_annot_04092016.tab.bgz | awk '\
  BEGIN {
    FS="\t"
    OFS="\t"
    print "##fileformat=VCFv4.2"
    print "##source=Eigen_hg19_coding_annot_04092016.tab.bgz"
    print "##INFO=<ID=EigenRaw,Number=1,Type=Float,Description=\"Eigen raw score\">"
    print "##INFO=<ID=EigenPhred,Number=1,Type=Float,Description=\"Eigen Phred-scaled score\">"
    print "##INFO=<ID=EigenPCRaw,Number=1,Type=Float,Description=\"EigenPC raw score\">"
    print "##INFO=<ID=EigenPCPhred,Number=1,Type=Float,Description=\"EigenPC Phred-scaled score\">"
    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  }

  (NR > 1) {
    print $1,$2,".",$3,$4,".",".","EigenRaw="$18";EigenPhred="$19";EigenPCRaw="$20";EigenPCPhred="$21
  }' | bgzip > ../converted_data/Eigen_hg19_coding_annot_04092016.vcf.gz

tabix ../converted_data/Eigen_hg19_coding_annot_04092016.vcf.gz

./bin/vt normalize -r ../../../resources/hs37d5x/reference_genome/hs37d5x.fa -o ../converted_data/Eigen_hg19_coding_annot_04092016.norm.vcf.gz ../converted_data/Eigen_hg19_coding_annot_04092016.vcf.gz

mv ../converted_data/Eigen_hg19_coding_annot_04092016.norm.vcf.gz ../converted_data/Eigen_hg19_coding_annot_04092016.norm.vcf.bgz
tabix ../converted_data/Eigen_hg19_coding_annot_04092016.norm.vcf.bgz
