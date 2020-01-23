#!/bin/bash
set -euo pipefail

# Command line arguments
GATK="$1"
REF_FA="$2"
LOCI_TSVXZ="$3"
INPUT_BAM="$4"
OUTPUT_AFSXZ="$5"

# Temporary files
loci_intervals=$(mktemp -s intervals)
genotype_vcf=$(mktemp -s vcf)

# Convert the TSV-format loci file into a GATK format .intervals
xz -dc "${LOCI_TSVXZ}" | awk '{print $1 ":" $2 "-" $2}' > "${loci_intervals}"

# Run GATK HC
# -ip 100 instructs HC to consider a region of 100 bp around each locus,
# to enable local haplotype reassembly.  Note that because of this, some 
# additional variant loci may be reported (not just those in LOCI_TSVXZ),
# but these will be removed at the later R stage.
java -Xmx2G -jar "${GATK}" -T HaplotypeCaller -R "${REF_FA}" -L "${loci_intervals}" -ip 100 -I "${INPUT_BAM}" -o "${genotype_vcf}"

# Post-process the GATK VCF: extract het SNP loci and report depths.
grep -v '^#' "${genotype_vcf}" | awk '
BEGIN {
    FS="\t"
    OFS="\t"
}

(
    ($4 == "G" || $4 == "C" || $4 == "A" || $4 == "T") && 
    ($5 == "G" || $5 == "C" || $5 == "A" || $5 == "T"))
{
    split($9, fields, ":")
    split($10, values, ":")
    for (i in fields)
        data[fields[i]] = values[i]

    if (data["GT"] != "0/1" && data["GT"] != "1/0")
        next

    split(data["AD"], ads, ",")
    print $1, $2, ads[0] + ads[1], ads[1]
}' | xz -c > "${OUTPUT_AFSXZ}"

# Clean up
rm -f "${loci_intervals}" "${genotype_vcf}"
