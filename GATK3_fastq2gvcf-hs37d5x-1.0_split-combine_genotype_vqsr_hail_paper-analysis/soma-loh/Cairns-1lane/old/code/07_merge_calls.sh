#!/bin/bash
set -euo pipefail

mkdir -p tmp

firstfile=0
for f in ../*.afs.aneu.tsv; do
  if [ "$firstfile" == "0" ]; then
    cp "$f" ./tmp/merged.txt
    firstfile=1
  else
    tail -n+2 "$f" >> ./tmp/merged.txt
  fi
done

mv ./tmp/merged.txt ../07_Cairns-1lane.afs.aneu.merged.tsv
