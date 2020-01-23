#!/usr/bin/env Rscript
options(error = stop, echo = TRUE)

suppressPackageStartupMessages(library(SomaticSignatures))
suppressPackageStartupMessages(library(plyr))

donefiles = list.files("../04_fact_runs/", pattern = "^04_combined_burden\\.fact\\.(single|agegrouped)\\.[0-9]+-[0-9]+-[0-9]+\\.[0-9]+-[0-9]+-[0-9]+\\.[0-9]+\\.rds\\.done$", full.names = TRUE)

descs = gsub(".*/04_combined_burden\\.fact\\.", "", gsub("\\.rds\\.done$", "", donefiles))
infiles = gsub("\\.done$", "", donefiles)
descparts = strsplit(descs, "[-.]")
descs = as.data.frame(t(simplify2array(descparts)), stringsAsFactors = FALSE)
colnames(descs) = c("mode", "kmin", "kmax", "B", "seed", "i", "N", "k")
for (i in 2:ncol(descs))
    descs[,i] = as.integer(descs[,i])
descs$infile = infiles
rm(descparts, donefiles, infiles)

fits_summary = ddply(descs, .(infile), function(d) {
        fit = readRDS(d$infile)
        d$evar = evar(fitted(fit), observed(fit))
        d$rss = rss(fitted(fit), observed(fit))
        d
    }, .progress = "text")

saveRDS(fits_summary, "../05_combined_factorise_summary.rds")


fits_best = dlply(fits_summary, .(mode, k), function(d) {
    readRDS(d$infile[which.max(d$evar)])
})

saveRDS(fits_best, "../05_combined_factorise_bestfits.rds")

