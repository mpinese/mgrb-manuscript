#!/bin/bash
set -euo pipefail

INFILE="../MGRB_phase2.dupmarked.realigned.recalibrated.cgc.freebayes.split.vtnorm.merged.tsv"

mkdir -p tmp

set +u
source bin/perlbrew_vep_environment.sh
set -u

# Generate an input file and feed to VEP
cut -f 2-5 "$INFILE" | \
  tail -n+2 | \
  LC_ALL=C sort -T /nvme/tmp/marpin --parallel=8 | \
  uniq | \
  awk 'BEGIN {OFS="\t"; print "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"} {print $1,$2,$1":"$2":"$3":"$4,$3,$4,".\t.\t."}' | \
  ./bin/vep \
    --format vcf -o tmp/vep.out --force --no_stats --warning_file tmp/vep.warnings \
    --assembly GRCh37 --minimal --gencode_basic \
    --variant_class --regulatory --total_length --numbers --no_escape \
    --hgvs --protein --symbol --canonical \
    --max_af \
    --tab \
    --no_intergenic --pick \
    --cache --offline --dir ../../software/.vep --fasta ../../software/.vep/homo_sapiens/90_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
    --fork 28

