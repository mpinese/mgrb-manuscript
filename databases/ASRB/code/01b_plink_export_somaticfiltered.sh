#!/bin/bash
set -euo pipefail

mkdir -p ../plink

./bin/pyhail-0.1-latest.sh ./01b_plink_export_somaticfiltered.py

cat ../asrb_male.sample_list ../asrb_female.sample_list | awk '{print 0, $1}' > tmp/samples.txt
bin/plink --bfile tmp/01b_genotypes \
  --keep tmp/samples.txt \
  --make-bed --freqx \
  --out ../plink/01b_ASRB.WGStier12.GiaB_HCR.split.minrep.somaticfiltered

awk '{print 0, $1}' < ../asrb_male.sample_list > tmp/samples.txt
bin/plink --bfile ../plink/01b_ASRB.WGStier12.GiaB_HCR.split.minrep.somaticfiltered \
  --keep tmp/samples.txt \
  --freqx \
  --out ../plink/01b_ASRB.WGStier12.GiaB_HCR.split.minrep.somaticfiltered.m

awk '{print 0, $1}' < ../asrb_female.sample_list > tmp/samples.txt
bin/plink --bfile ../plink/01b_ASRB.WGStier12.GiaB_HCR.split.minrep.somaticfiltered \
  --keep tmp/samples.txt \
  --freqx \
  --out ../plink/01b_ASRB.WGStier12.GiaB_HCR.split.minrep.somaticfiltered.f


