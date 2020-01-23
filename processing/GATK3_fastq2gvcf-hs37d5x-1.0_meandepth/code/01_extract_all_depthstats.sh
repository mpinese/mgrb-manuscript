#!/bin/bash
for f in ../../GATK3_fastq2gvcf-hs37d5x-1.0/bams/misc/*.misc.tar.gz; do
  echo $f
  ./extract_depthstats.sh "$f" ..
done

for f in ../*.stats; do echo -e $(basename $f | sed 's/_.*//')\\t$(grep MEAN $f | sed 's/.* //'); done > ../meandepth.txt
