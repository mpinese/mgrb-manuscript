#!/bin/bash
set -euo pipefail

mkdir -p ../plink

./bin/pyhail-0.1-latest.sh ./08b_plink_export_somaticfiltered.py

rm -f tmp/bedlist.txt
seq 2 22 | while read i; do
  echo "tmp/08b_genotypes_${i}.bed tmp/08b_genotypes_${i}.bim tmp/08b_genotypes_${i}.fam" >> tmp/bedlist.txt
done

bin/plink --bfile tmp/08b_genotypes_1 --merge-list tmp/bedlist.txt \
  --make-bed --freqx \
  --out ../plink/08b_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.GiaB_HCR.split.minrep.somaticfiltered

awk '{print 0, $1}' < ../mgrb_phase2_male.sample_list > tmp/samples.txt
bin/plink --bfile ../plink/08b_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.GiaB_HCR.split.minrep.somaticfiltered \
  --keep tmp/samples.txt \
  --freqx \
  --out ../plink/08b_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.m.GiaB_HCR.split.minrep.somaticfiltered

awk '{print 0, $1}' < ../mgrb_phase2_female.sample_list > tmp/samples.txt
bin/plink --bfile ../plink/08b_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.GiaB_HCR.split.minrep.somaticfiltered \
  --keep tmp/samples.txt \
  --freqx \
  --out ../plink/08b_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.f.GiaB_HCR.split.minrep.somaticfiltered

awk '{print 0, $1}' < ../mgrb_phase2_noch10.sample_list > tmp/samples.txt
bin/plink --bfile ../plink/08b_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.GiaB_HCR.split.minrep.somaticfiltered \
  --keep tmp/samples.txt \
  --freqx \
  --out ../plink/08b_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.noch10.GiaB_HCR.split.minrep.somaticfiltered

awk '{print 0, $1}' < ../mgrb_phase2_ch10.sample_list > tmp/samples.txt
bin/plink --bfile ../plink/08b_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.GiaB_HCR.split.minrep.somaticfiltered \
  --keep tmp/samples.txt \
  --freqx \
  --out ../plink/08b_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.ch10.GiaB_HCR.split.minrep.somaticfiltered


