#!/bin/bash
curl http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz | \
  gzip -dc | \
  cut -f 1-3 | \
  sed 's/^chr//; s/^M/MT/' | \
  sort -k1,1 -k2,2n \
    > wgEncodeDacMapabilityConsensusExcludable.bed
