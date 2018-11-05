#!/bin/bash
set -euo pipefail

plink \
  --bfile ../plink/08_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.GiaB_HCR.split.minrep \
  --remove 09_dreamlab_dropped_samples.txt \
  --geno 0.01 \
  --maf 0.01 \
  --hwe 1e-6 midp \
  --pca header tabs --make-bed \
  --out ../plink/09_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.GiaB_HCR.split.minrep.fordreamlab
