#!/bin/bash
# Note: Run
# perlbrew use perl-5.26.0
# in shell before executing

# Prepare input file for VEP
mkdir -p temp
awk 'BEGIN{FS=",";OFS="\t";print "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"} {print $1, $2, $1 "," $2 "," $3 "," $4, $3, $4, ".\t.\t."}' < ../../allelefreqs/MGRB_GnomAD_SweGen_cancer_AFs_sharedhcregions_ccdsexonsplus2bp.csv > temp/vep_input.vcf

set -euo pipefail

../../../software/ensembl-vep/vep \
  --format vcf -i temp/vep_input.vcf \
  --tab -o temp/vep_output.tab --no_stats --force_overwrite \
  --assembly GRCh37 --minimal --everything --pick --coding_only --gencode_basic \
  --cache --offline --dir ../../../software/.vep --fasta ../../../software/.vep/homo_sapiens/90_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
  --fork 28

# Filter VEP results: keep only hits in canonical protein coding transcripts with HGNC IDs, and return a subset of fields
grep -v '^##' temp/vep_output.tab \
  | sed 's/^#//' \
  | awk 'BEGIN{FS="\t";OFS="\t"} (($25 == "YES" && $22 == "HGNC" && $24 == "protein_coding") || NR == 1) { if ($28 == "-") { $28 = "" }; print $1,$5,$7,$10,$18,$21,$28,$36,$40 }' \
  | sed 's/%3D/=/g; s/splice_region_variant//; s/start_retained_variant//; s/stop_retained_variant//; s/3_prime_UTR_variant//; s/5_prime_UTR_variant//; s/coding_sequence_variant//' \
  | sed 's/,+/,/g; s/\t,/\t/; s/,\t/\t/' \
  | awk 'BEGIN{FS="\t";OFS="\t"} (NR==1) { $1="contig\tpos\tref\talt";print } (NR > 1 && $3 != "") { gsub(/,/, "\t", $1); print}' \
  > ../MGRB_GnomAD_SweGen_cancer_AFs_sharedhcregions_ccdsexonsplus2bp.vep.filtered.tsv



#   ColumnIndex    Desc
#   1              Uploaded_variation : Identifier of uploaded variant
#   5              Feature : Stable ID of feature
#   6              Feature_type : Type of feature - Transcript, RegulatoryFeature or MotifFeature
#   7              Consequence : Consequence type
#   10             Protein_position : Relative position of amino acid in protein
#   18             FLAGS : Transcript quality flags
#   21             SYMBOL : Gene symbol (e.g. HGNC)
#   22             SYMBOL_SOURCE : Source of gene symbol
#   23             HGNC_ID : Stable identifer of HGNC gene symbol
#   24             BIOTYPE : Biotype of transcript or regulatory feature
#   25             CANONICAL : Indicates if transcript is canonical for this gene
#   28             CCDS : Indicates if transcript is a CCDS transcript
#   36             EXON : Exon number(s) / total
#   40             HGVSp : HGVS protein sequence name

