#!/usr/bin/env Rscript
options(error = stop)

suppressPackageStartupMessages(library(SomaticSignatures))
suppressPackageStartupMessages(library(plyr))

donefiles = list.files("..", pattern = "^02_lowdepth_burden.fact.[0-9]+-[0-9]+-[0-9]+\\.[0-9]+-[0-9]+-[0-9]+\\.[0-9]+\\.rds\\.done$", full.names = TRUE)

descs = gsub(".*/02_lowdepth_burden\\.fact\\.", "", gsub("\\.rds\\.done$", "", donefiles))
infiles = gsub("\\.done$", "", donefiles)
descparts = strsplit(descs, "[-.]")
descs = as.data.frame(t(simplify2array(descparts)), stringsAsFactors = FALSE)
colnames(descs) = c("kmin", "kmax", "B", "seed", "i", "N", "k")
for (i in 1:ncol(descs))
    descs[,i] = as.integer(descs[,i])
descs$infile = infiles
rm(descparts, donefiles, infiles)

fits_summary = ddply(descs, .(infile), function(d) {
        fit = readRDS(d$infile)
        d$evar = evar(fitted(fit), observed(fit))
        d$rss = rss(fitted(fit), observed(fit))
        d
    }, .progress = "text")

saveRDS(fits_summary, "../03_lowdepth_factorise_summary.rds")


fits_best = dlply(fits_summary, .(k), function(d) {
    readRDS(d$infile[which.max(d$evar)])
})

saveRDS(fits_best, "../03_lowdepth_factorise_bestfits.rds")

