#!/usr/bin/env Rscript

pdf("../02_calib.pdf", height = 8, width = 8)	

message("Reading MGRB depths")
dp.mgrb = read.table("../../MGRB/03_dp_stats.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
message("Reading Cairns depths")
dp.cairns = read.table("../../Cairns-1lane/05_dp_stats.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
message("Reading GC")
gc = read.table("../unfiltered/commonhets.gc.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
message("Reading variants")
vars = read.table("../unfiltered/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.commonhets.variant_list", header = FALSE, sep = ":", stringsAsFactors = FALSE)
colnames(vars) = c("chrom", "pos", "ref", "alt")

colnames(dp.mgrb) = c("chrom", "pos", "meandp.mgrb", "vardp.mgrb")
colnames(dp.cairns) = c("chrom", "pos", "meandp.cairns", "vardp.cairns")

message("Merging")
merged = merge(merge(merge(vars, gc), dp.mgrb), dp.cairns)

rm(dp.mgrb, dp.cairns, gc, vars)

message("Filtering")
pairs(merged[1:10000,c("gc100", "gc400", "meandp.mgrb", "meandp.cairns", "vardp.mgrb", "vardp.cairns")], pch = ".")

merged = merged[merged$gc100 >= 0.3 & merged$gc100 <= 0.55,]

k.mgrb = 1/mean(merged$meandp.mgrb)
k.cairns = 1/mean(merged$meandp.cairns)
merged$meandp.mgrb = merged$meandp.mgrb * k.mgrb
merged$meandp.cairns = merged$meandp.cairns * k.cairns
merged$vardp.mgrb = merged$vardp.mgrb * k.mgrb^2
merged$vardp.cairns = merged$vardp.cairns * k.cairns^2

pairs(merged[1:10000,c("gc100", "gc400", "meandp.mgrb", "meandp.cairns", "vardp.mgrb", "vardp.cairns")], pch = ".")

merged = merged[
    merged$vardp.mgrb >= 0.025 & merged$vardp.mgrb <= 0.033 & 
    merged$vardp.cairns >= 0.025 & merged$vardp.cairns < 0.04 & 
    merged$meandp.mgrb >= 0.9 & merged$meandp.mgrb <= 1.1 &
    merged$meandp.cairns >= 0.9 & merged$meandp.cairns <= 1.1,]

k.mgrb = 1/mean(merged$meandp.mgrb)
k.cairns = 1/mean(merged$meandp.cairns)
merged$meandp.mgrb = merged$meandp.mgrb * k.mgrb
merged$meandp.cairns = merged$meandp.cairns * k.cairns
merged$vardp.mgrb = merged$vardp.mgrb * k.mgrb^2
merged$vardp.cairns = merged$vardp.cairns * k.cairns^2

merged$meandp.delta = (merged$meandp.mgrb - merged$meandp.cairns)

pairs(merged[1:10000,c("gc100", "gc400", "meandp.mgrb", "meandp.cairns", "vardp.mgrb", "vardp.cairns", "meandp.delta")], pch = ".")

merged = merged[abs(merged$meandp.delta) < 0.05,]

k.mgrb = 1/mean(merged$meandp.mgrb)
k.cairns = 1/mean(merged$meandp.cairns)
merged$meandp.mgrb = merged$meandp.mgrb * k.mgrb
merged$meandp.cairns = merged$meandp.cairns * k.cairns
merged$vardp.mgrb = merged$vardp.mgrb * k.mgrb^2
merged$vardp.cairns = merged$vardp.cairns * k.cairns^2

pairs(merged[1:10000,c("gc100", "gc400", "meandp.mgrb", "meandp.cairns", "vardp.mgrb", "vardp.cairns", "meandp.delta")], pch = ".")

message("Writing output")
write.table(merged[,c("chrom", "pos", "ref", "alt")], file = "../soma-loh.hs37.loci.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(merged[,c("chrom", "pos", "gc100", "gc200", "gc400", "gc600", "gc800")], file = "../soma-loh.hs37.gc.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

temp = merged[,c("chrom", "pos", "meandp.mgrb")]
colnames(temp)[3] = "affinity"
write.table(temp, file = "../soma-loh.hs37.affinity.mgrb.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

temp = merged[,c("chrom", "pos", "meandp.cairns")]
colnames(temp)[3] = "affinity"
write.table(temp, file = "../soma-loh.hs37.affinity.cairns.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

dev.off()

