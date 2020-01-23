#!/bin/bash
set -euo pipefail

find .. -maxdepth 1 -type d -name '[ABZ][A-Z][A-Z][A-Z][A-Z]' | ./bin/parallel --gnu Rscript ./03h_load_variants.onesample.R {}

Rscript ./03h_load_variants.combinesamples.R
