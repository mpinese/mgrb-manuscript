#!/bin/bash
set -euo pipefail

first_file="1"
find .. -maxdepth 1 -name '*.cgc.freebayes.split.vtnorm.tsv' | sort | while read path; do
  echo "${path}" > /dev/stderr
  if [ "${first_file}" == "1" ]; then
    head -n1 "${path}" > "../MGRB_phase2.dupmarked.realigned.recalibrated.cgc.freebayes.split.vtnorm.merged.tsv" || true
    first_file="0"
  fi

  tail -n+2 "${path}" >> "../MGRB_phase2.dupmarked.realigned.recalibrated.cgc.freebayes.split.vtnorm.merged.tsv"
done
