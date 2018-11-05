#!/bin/bash
set -euo pipefail

# Usage: extract_depthstats.sh <targzfile> <destdir>
sampleid=$(basename "$1")
sampleid="${sampleid:0:5}"

ds_file=$(tar tzf "$1" | grep "${sampleid}_.*\\.dupmarked.realigned.recalibrated.bam.depth.stats")

tar xzf "$1" --strip-components=3 -C "$2" "$ds_file"
