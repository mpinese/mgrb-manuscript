#!/bin/bash
set -euo pipefail

mkdir -p ../source_data
wget https://raw.githubusercontent.com/joerivandervelde/gavin/master/data/predictions/GAVIN_ruleguide_r0.3.tsv -O ../source_data/GAVIN_ruleguide_r0.3.tsv
echo $(date) >> ../source_data/GAVIN_ruleguide_r0.3.download_timestamp
