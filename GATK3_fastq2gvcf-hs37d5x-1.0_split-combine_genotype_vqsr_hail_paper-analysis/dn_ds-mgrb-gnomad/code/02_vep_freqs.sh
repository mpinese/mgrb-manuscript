#!/bin/bash

# Prepare input files for VEP
gzip -dc ../01_gnomad.genomes.r2.0.1.sites.autosomes.split.all.freqs.vcf.bgz | \
  grep -vF '#' | grep -vF 'AC_NFE=0;' | sed -E 's/;/\t/g; s/[A-Za-z]+_NFE=//g' | \
  awk 'BEGIN{FS="\t";OFS="\t";print "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tINFO"} {print $1, $2, $1 ":" $2 ":" $4 ":" $5 "_GnomAD:" $9 ":" $10 ":" $11, $4, $5, ".\t.\t."}' | \
  bgzip > ../02_gnomad.genomes.r2.0.1.sites.autosomes.split.all.freqs.minimal.vcf.gz

gzip -dc ../01_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.all.freqs.vcf.bgz | \
  grep -vF '#' | grep -vF 'AC_MGRB=0;' | sed -E 's/;/\t/g; s/[A-Za-z]+_MGRB=//g' | \
  awk 'BEGIN{FS="\t";OFS="\t";print "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"} {print $1, $2, $1 ":" $2 ":" $4 ":" $5 "_MGRB:" $9 ":" $10 ":" $11, $4, $5, ".\t.\t."}' | \
  bgzip > ../02_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.all.freqs.minimal.vcf.gz

# Run VEP
bash ./02h_vep_helper.sh

# Filter VEP results: keep only hits in canonical protein coding transcripts with HGNC IDs, and return a subset of fields
grep -v '^##' ../02_gnomad.genomes.r2.0.1.sites.autosomes.split.all.freqs.minimal.vep.tab | sed 's/^#//' | \
  awk 'BEGIN{FS="\t";OFS="\t"} (($25 == "YES" && $22 == "HGNC" && $24 == "protein_coding") || NR == 1) { if ($18 == "-") { $18 = "" }; if ($28 == "-") { $28 = "" }; print $1,$5,$7,$18,$21,$28,$36,$40 }' > \
  ../02_gnomad.genomes.r2.0.1.sites.autosomes.split.all.freqs.minimal.vep.filtered.tsv

grep -v '^##' ../02_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.all.freqs.minimal.vep.tab | sed 's/^#//' | \
  awk 'BEGIN{FS="\t";OFS="\t"} (($25 == "YES" && $22 == "HGNC" && $24 == "protein_coding") || NR == 1) { if ($18 == "-") { $18 = "" }; if ($28 == "-") { $28 = "" }; print $1,$5,$7,$18,$21,$28,$36,$40 }' > \
  ../02_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.all.freqs.minimal.vep.filtered.tsv

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

