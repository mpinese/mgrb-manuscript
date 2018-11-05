#!/bin/bash
set -euo pipefail

mkdir -p ../plink
mkdir -p ./tmp

./bin/pyhail.sh ./01_plink_export.py

cat ../asrb_male.sample_list ../asrb_female.sample_list | awk '{print 0, $1}' > tmp/samples.txt
bin/plink --bfile tmp/01_genotypes \
  --keep tmp/samples.txt \
  --make-bed --freqx \
  --out ../plink/01_ASRB.WGStier12.GiaB_HCR.split.minrep


awk '{print 0, $1}' < ../asrb_male.sample_list > tmp/samples.txt
bin/plink --bfile ../plink/01_ASRB.WGStier12.GiaB_HCR.split.minrep \
  --keep tmp/samples.txt \
  --freqx \
  --out ../plink/01_ASRB.WGStier12.GiaB_HCR.split.minrep.m


awk '{print 0, $1}' < ../asrb_female.sample_list > tmp/samples.txt
bin/plink --bfile ../plink/01_ASRB.WGStier12.GiaB_HCR.split.minrep \
  --keep tmp/samples.txt \
  --freqx \
  --out ../plink/01_ASRB.WGStier12.GiaB_HCR.split.minrep.f
