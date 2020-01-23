#!/usr/bin/env Rscript

# Load data
temp.cairns = readRDS("../input/low-depth/03_cairns_soma-snv_pass2.combined.burden.rds")
temp.mgrb = readRDS("../input/low-depth/03_mgrb_soma-snv_pass2.combined.burden.rds")


# Keep only high quality samples
temp.mgrb.tier1 = read.table("../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/tier1.sample_list", stringsAsFactors = FALSE)[,1]
temp.mgrb.tier2 = read.table("../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/tier2.sample_list", stringsAsFactors = FALSE)[,1]
temp.mgrb$count = temp.mgrb$count[,colnames(temp.mgrb$count) %in% c(temp.mgrb.tier1, temp.mgrb.tier2)]
temp.mgrb$background = temp.mgrb$background[,colnames(temp.mgrb$background) %in% c(temp.mgrb.tier1, temp.mgrb.tier2)]
stopifnot(rownames(temp.mgrb$count) == rownames(temp.mgrb$background))
stopifnot(colnames(temp.mgrb$count) == colnames(temp.mgrb$background))

temp.cairns.good = read.table("../../soma-loh/Cairns-1lane/04_good_samples.txt", stringsAsFactors = FALSE)[,1]
# Cairns-1-FR07888505-H75N5CCXX-7 is a conspicious outlier, with ten times
# the burden of the next highest burden sample.  Remove it for now for the
# sake of the fits, and later revisit to determine the cause.
temp.cairns.good = setdiff(temp.cairns.good, "Cairns-1-FR07888505-H75N5CCXX-7")
temp.cairns$count = temp.cairns$count[,colnames(temp.cairns$count) %in% temp.cairns.good]
temp.cairns$background = temp.cairns$background[,colnames(temp.cairns$background) %in% temp.cairns.good]
stopifnot(rownames(temp.cairns$count) == rownames(temp.cairns$background))
stopifnot(colnames(temp.cairns$count) == colnames(temp.cairns$background))


# Merge burden matrices and rescale to mutations per megabase
stopifnot(rownames(temp.mgrb$count) == rownames(temp.cairns$count))
counts = cbind(temp.mgrb$count, temp.cairns$count)
background = cbind(temp.mgrb$background, temp.cairns$background)
burden = counts / background * 1e6
rm(list = grep("^temp", ls(), value = TRUE))


# Write the burden matrix
args = commandArgs(trailingOnly = TRUE)
saveRDS(burden, args[[1]])
