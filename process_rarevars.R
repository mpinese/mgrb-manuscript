# library(openxlsx)
# tbl = read.xlsx("data/Mark_and_David_variants_for_MGRB_paper_annotated_all_including_NANS_TRIM14_spreadsheet_extra_fields.xlsx")
# tbl = read.xlsx("data/Mark_and_David_variants_for_MGRB_paper_annotated_all_spreadsheet_extra_fields_FIXED_CONDEL_COLUMNS.xlsx")
tbl = read.table("data/Mark_David_for_MGRB_paper_2018oct_all_spreadsheet_extra_fields.tsv.xz", header = TRUE, stringsAsFactors = FALSE, sep = "\t")


# Initial filtering (File is already SNVs only):
# Coding or splice site only
# AF < 5% in any cohort in MGRB, GnomAD NFE, SweGen
# In comparison regions

# Filter by coding variants only
# Note we include synonymous vars here.
# We exclude some variants which are of interest, but
# likely under very different constraints:
# splice_acceptor_variant splice_donor_variant start_lost stop_lost
tbl = tbl[tbl$vep_consequence_terms %in% c("missense_variant", "stop_gained", "synonymous_variant"),]

# Drop unneeded or uninformative columns
tbl = tbl[,!(colnames(tbl) %in% c("VAF", "cohort_MAF", "is_tier3variant", "has_POT1_as_interactor_A", "has_POT1_as_interactor_B", "gnomad_AF_NFE", "gnomad_AF", "swegen_AF"))]

# Parse variant information from the VID
temp = strsplit(tbl$VARIANT, ":")
tbl$contig = as.factor(sapply(temp, function(x) x[1]))
tbl$pos = as.integer(sapply(temp, function(x) x[2]))
tbl$ref = as.factor(sapply(temp, function(x) x[3]))
tbl$alt = as.factor(sapply(temp, function(x) x[4]))

# Subset to SNVs
tbl = tbl[tbl$ref %in% c("A", "C", "G", "T") & tbl$alt %in% c("A", "C", "G", "T"),]

# Subset to variants in the comparison regions
library(GenomicRanges)
temp = read.table("data/MGRBphase2final_dpge15_98pct.GnomAD201_dpge15_98pct.GiaB332HC.ccdsexonsplus2bp.bed.xz", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
# 1-based inclusive:
comparison_regions = GRanges(
    seqnames = Rle(as.character(temp[,1])),
    ranges = IRanges(start = temp[,2]+1, end = temp[,3]),
    strand = rep('*', nrow(temp)))
variant_regions = GRanges(
    seqnames = Rle(as.character(tbl$contig)),
    ranges = IRanges(start = tbl$pos, width = 1),
    strand = rep('*', nrow(tbl)),
    vid = tbl$VARIANT)
variant_comparison_region_intersect = countOverlaps(variant_regions, comparison_regions)
tbl = tbl[tbl$VARIANT %in% variant_regions$vid[variant_comparison_region_intersect != 0],]


# Convert columns to aid in compression
for (i in 1:ncol(tbl))
{
    if (class(tbl[,i]) == "character" && mean(table(tbl[,i])) > 1.7)
    {
        message(sprintf("Converting %s to factor", colnames(tbl)[i]))
        tbl[,i] = factor(tbl[,i])
    }
    else if (grepl("^n(RR|RA|AA|missing)_.*", colnames(tbl)[i]) && class(tbl[,i]) == "numeric")
    {
        message(sprintf("Converting %s to integer", colnames(tbl)[i]))
        tbl[,i] = as.integer(tbl[,i])
    }
}
tbl$CosmicCodingMuts_CNT = as.integer(tbl$CosmicCodingMuts_CNT)
tbl$CosmicMutantExport_HGNCId = as.integer(tbl$CosmicMutantExport_HGNCId)


# Add RVIS scores for genes
# RVIS scores are from genic-intolerance.org, location
# http://genic-intolerance.org/data/RVIS_Unpublished_ExACv2_March2017.txt
# description:
# "novel unpublished RVIS gene score based on ExAC v2 release 2.0 (accessed: March 15th 2017). 
# As of this release we use CCDS release 20 and Ensembl release 87 annotations."
# Downloaded 5 Oct 2018.
rvis = read.table("data/RVIS_Unpublished_ExACv2_March2017.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
temp.rvis_match = match(as.character(tbl$gene_symbol), rvis$CCDSr20)
tbl$rvis.coverage = rvis$X.geneCov[temp.rvis_match]
tbl$rvis.rvis_percentile = rvis$X.RVIS.pop_maf_0.05..any..[temp.rvis_match]
tbl$rvis.edge_case = rvis$Edge_case_RVIS.pop_maf_0.05..any..[temp.rvis_match] == "Y"


# Add old (pre gnomAD genomes public release) ClinVar scores
clinvar = read.table("data/variant_summary_2016-09.txt.gz", sep = "\t", header = FALSE, quote = "", skipNul = TRUE, stringsAsFactors = FALSE)
# SNVs on hg37 only
clinvar = clinvar[clinvar$V13 == "GRCh37" & clinvar$V15 == clinvar$V16 & clinvar$V26 %in% c("G", "C", "A", "T") & clinvar$V27 %in% c("G", "C", "A", "T"),]
clinvar$VARIANT = sprintf("%s:%d:%s:%s", clinvar$V14, clinvar$V15, clinvar$V26, clinvar$V27)
clinvar = clinvar[,c("VARIANT", "V6", "V18")]
colnames(clinvar) = c("VARIANT", "Clinvar201609.sig", "Clinvar201609.evidence")
clinvar = clinvar[!duplicated(clinvar),]
tbl = merge(tbl, clinvar, by = "VARIANT", all.x = TRUE, all.y = FALSE)

tbl[tbl == ""] = NA

saveRDS(tbl, "data/rare_coding_varfreqs_mgrb_gnomad_swegen.rds")
