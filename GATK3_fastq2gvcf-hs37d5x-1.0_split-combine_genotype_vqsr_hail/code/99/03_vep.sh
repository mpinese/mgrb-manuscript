#!/bin/bash

parallel ../../software/ensembl-vep/vep -v -i {} --format vcf -o {}.vep.txt --tab --everything --hgvsg --flag_pick --no-escape --database ../../software/.vep/ --fasta ../../../resources/hs37d5x/reference_genome/hs37d5x.fa ::: MGRB.phase2.tier12.match.vqsr.minrep.split.variants.vcf.bgz/*.bgz

../../software/ensembl-vep/vep -v -i {} --format  -o ../acmg.vep.vcf --vcf --everything --hgvsg --flag_pick_allele --no_escape --offline --dir_cache ../../software/.vep/ --fasta ../../../resources/hs37d5x/reference_genome/hs37d5x.fa
