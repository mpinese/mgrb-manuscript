#!/usr/bin/env Rscript

infiles = list.files("..", ".*\\.afs\\.bgz$", full.names = TRUE)

library(parallel)

depths = data.frame(t(simplify2array(mclapply(infiles, function(infile) fivenum(scan(infile, what = list(NULL, NULL, integer(), NULL))[[3]]), mc.cores = detectCores()))))
depths = cbind(gsub("\\.afs\\.bgz$", "", gsub(".*/", "", infiles)), depths)
colnames(depths) = c("sample", "min", "h1", "med", "h2", "max")
depths$sample = as.character(depths$sample)

write.table(depths, "../03_Cairns-1lane.median_depths.tsv", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

