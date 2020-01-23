#!/bin/bash
set -euo pipefail

mkdir -p ../plink

./bin/pyhail.sh ./07_plink_export_45andup_subpops.py

for f in ../plink/07*.bed; do
  ./bin/plink --bfile "${f%.bed}" --freqx --out "${f%.bed}"
done

