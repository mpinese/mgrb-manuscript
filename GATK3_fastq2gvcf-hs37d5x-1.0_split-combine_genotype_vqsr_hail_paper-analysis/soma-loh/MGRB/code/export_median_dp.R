#!/usr/bin/env Rscript

infiles = list.files("../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.commonhets", ".*\\.afs\\.bgz$", full.names = TRUE)

library(parallel)

depths = data.frame(t(simplify2array(mclapply(infiles, function(infile) fivenum(scan(infile, what = list(NULL, NULL, integer(), NULL))[[3]]), mc.cores = detectCores()))))
depths = cbind(gsub("\\.afs\\.bgz$", "", gsub(".*/", "", infiles)), depths)
colnames(depths) = c("sample", "min", "h1", "med", "h2", "max")
depths$sample = as.character(depths$sample)

write.table(depths, "../MGRB.median_depths.tsv", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

