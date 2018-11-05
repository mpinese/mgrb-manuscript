#!/usr/bin/env Rscript
#options(echo = TRUE)

# LOAD LOW DEPTH DATA
############################################

lowdepth.cairns = readRDS("../input/low-depth/03_cairns_soma-snv_pass2.combined.burden.rds")
lowdepth.mgrb = readRDS("../input/low-depth/03_mgrb_soma-snv_pass2.combined.burden.rds")
lowdepth.swegen = readRDS("../input/low-depth/03_swegen_soma-snv_pass2.combined.burden.rds")

# Keep only high quality samples
lowdepth.mgrb.tier1 = read.table("../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/tier1.sample_list", stringsAsFactors = FALSE)[,1]
lowdepth.mgrb.tier2 = read.table("../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/tier2.sample_list", stringsAsFactors = FALSE)[,1]
lowdepth.mgrb$count = lowdepth.mgrb$count[,colnames(lowdepth.mgrb$count) %in% c(lowdepth.mgrb.tier1, lowdepth.mgrb.tier2)]
lowdepth.mgrb$background = lowdepth.mgrb$background[,colnames(lowdepth.mgrb$background) %in% c(lowdepth.mgrb.tier1, lowdepth.mgrb.tier2)]
stopifnot(rownames(lowdepth.mgrb$count) == rownames(lowdepth.mgrb$background))
stopifnot(colnames(lowdepth.mgrb$count) == colnames(lowdepth.mgrb$background))

lowdepth.cairns.good = read.table("../../soma-loh/Cairns-1lane/04_good_samples.txt", stringsAsFactors = FALSE)[,1]
# Cairns-1-FR07888505-H75N5CCXX-7 is a conspicious outlier, with ten times
# the burden of the next highest burden sample.  Remove it for now for the
# sake of the fits, and later revisit to determine the cause.
lowdepth.cairns.good = setdiff(lowdepth.cairns.good, "Cairns-1-FR07888505-H75N5CCXX-7")
lowdepth.cairns$count = lowdepth.cairns$count[,colnames(lowdepth.cairns$count) %in% lowdepth.cairns.good]
lowdepth.cairns$background = lowdepth.cairns$background[,colnames(lowdepth.cairns$background) %in% lowdepth.cairns.good]
stopifnot(rownames(lowdepth.cairns$count) == rownames(lowdepth.cairns$background))
stopifnot(colnames(lowdepth.cairns$count) == colnames(lowdepth.cairns$background))

# SweGen0315 is suspicious with 36 x the signal of the next highest
lowdepth.swegen$count = lowdepth.swegen$count[,colnames(lowdepth.swegen$count) != "SweGen0315"]
lowdepth.swegen$background = lowdepth.swegen$background[,colnames(lowdepth.swegen$background) != "SweGen0315"]

# Add LD- suffix to sample IDs
colnames(lowdepth.mgrb$count) = paste("MGRB-", colnames(lowdepth.mgrb$count), "-LD", sep = "")
colnames(lowdepth.mgrb$background) = paste("MGRB-", colnames(lowdepth.mgrb$background), "-LD", sep = "")
colnames(lowdepth.cairns$count) = paste(gsub("-[A-Z0-9]{9}-[0-8]$", "", gsub("^Cairns-1-", "Cairns-", colnames(lowdepth.cairns$count))), "-LD", sep = "")
colnames(lowdepth.cairns$background) = paste(gsub("-[A-Z0-9]{9}-[0-8]$", "", gsub("^Cairns-1-", "Cairns-", colnames(lowdepth.cairns$background))), "-LD", sep = "")
colnames(lowdepth.swegen$count) = paste("SweGen-", gsub("SweGen", "", colnames(lowdepth.swegen$count)), "-LD", sep = "")
colnames(lowdepth.swegen$background) = paste("SweGen-", gsub("SweGen", "", colnames(lowdepth.swegen$background)), "-LD", sep = "")


# LOAD HIGH DEPTH DATA
############################################

highdepth.all = readRDS("../input/high-depth/03_highdepth_soma-snv_pass2.combined.burden.rds")

# Dubbo 1002_1 is a known bad (admixed?) sample
highdepth.all$count = highdepth.all$count[,colnames(highdepth.all$count) != "Dubbo-1002_1"]
highdepth.all$background = highdepth.all$background[,colnames(highdepth.all$background) != "Dubbo-1002_1"]

# Add HD- prefix to sample IDs
colnames(highdepth.all$count) = paste(gsub("^ASPREE-", "MGRB-", colnames(highdepth.all$count)), "-HD", sep = "")
colnames(highdepth.all$background) = paste(gsub("^ASPREE-", "MGRB-", colnames(highdepth.all$background)), "-HD", sep = "")



# MERGE AND RESCALE TO MUTS / MB
############################################
stopifnot(rownames(lowdepth.mgrb$count) == rownames(lowdepth.cairns$count))
stopifnot(rownames(lowdepth.mgrb$count) == rownames(lowdepth.swegen$count))
stopifnot(rownames(lowdepth.mgrb$count) == rownames(highdepth.all$count))
stopifnot(rownames(lowdepth.mgrb$background) == rownames(lowdepth.cairns$background))
stopifnot(rownames(lowdepth.mgrb$background) == rownames(lowdepth.swegen$background))
stopifnot(rownames(lowdepth.mgrb$background) == rownames(highdepth.all$background))

counts = cbind(lowdepth.mgrb$count, lowdepth.cairns$count, lowdepth.swegen$count, highdepth.all$count)
background = cbind(lowdepth.mgrb$background, lowdepth.cairns$background, lowdepth.swegen$background, highdepth.all$background)
burden = counts / background * 1e6

# Write the burden matrix
args = commandArgs(trailingOnly = TRUE)
saveRDS(burden, args[[1]])
saveRDS(burden[,grepl("-LD$", colnames(burden))], args[[2]])




# SUMMARY STATISTICS
############################################
library(reshape2)
library(plyr)
burden = readRDS("../04_combined_burden.single.hdld.rds")
burden.melted = melt(burden, varnames = c("motif", "sample"), value.name = "burden")
burden.melted$cohort = gsub("-.*", "", burden.melted$sample)
burden.melted$depth = gsub(".*-", "", burden.melted$sample)
ddply(burden.melted, .(cohort, depth), function(d) fivenum(d$burden))
#   cohort depth V1 V2         V3         V4         V5
# 1 Cairns    HD  0  0 0.00000000 0.05527286  17.718560
# 2 Cairns    LD  0  0 0.00000000 0.00000000  66.342729
# 3  Dubbo    HD  0  0 0.03027797 0.09742165   9.188604
# 4   MGRB    HD  0  0 0.04235869 0.14715111  11.064014
# 5   MGRB    LD  0  0 0.00000000 0.05385006  26.258090
# 6 SweGen    LD  0  0 0.09183340 0.31787511 123.819435




# COLLAPSE BY AGE GROUP AND COHORT
############################################
# To reduce noise, produce a collapsed burden matrix, in which
# groups of ~16 samples of similar age from one cohort are combined.
# Do this only for the LD samples, due to the limited numbers of HD.
burden.ld = burden[,!grepl("-HD$", colnames(burden))]
# Get ages
ages = read.csv("ages.csv", stringsAsFactors = FALSE)
burden.ld.age = ages$AgeAtCollectionYears[match(colnames(burden.ld), ages$sampleID)]
# Group samples into sets of approximately 16 by cohort and age.
grouping = data.frame(sample = colnames(burden.ld), cohort = gsub("-.*", "", colnames(burden.ld)), age = burden.ld.age)
library(plyr)
grouping = ddply(grouping, .(cohort), function(d) {
    d = d[order(d$age),]
    numgroups = round(nrow(d)/16)
    d$group = paste(d$cohort, as.integer(floor(seq(1, numgroups + 1, length.out = nrow(d) + 1)[1:nrow(d)])), sep = "-G")
    d})

burden.ld.collapsed = t(daply(grouping, .(group), function(d) rowSums(burden.ld[,d$sample])))

#sort(colSums(burden.ld.collapsed))
burden.ld.collapsed.normalized = t(t(burden.ld.collapsed) / colSums(burden.ld.collapsed))

# Write the collapsed and normalized burden matrix
saveRDS(burden.ld.collapsed.normalized, args[[3]])
