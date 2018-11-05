#!/usr/bin/env Rscript
library(doParallel)
suppressPackageStartupMessages(library(SomaticSignatures))

# Fetch environment specs
env = Sys.getenv(c("INRDS", "OUTPREFIX", "SEED", "KMIN", "KMAX", "B", "CORES"))

seed = as.integer(env[["SEED"]])
kmin = as.integer(env[["KMIN"]])
kmax = as.integer(env[["KMAX"]])
B = as.integer(env[["B"]])
cores = as.integer(env[["CORES"]])

burden = readRDS(env[["INRDS"]])

zeroes = apply(burden == 0, 2, all)
burden = burden[,!zeroes]

# Factorise
for (i in kmin:kmax) {
    outfile = sprintf("%s.%d.rds", env[["OUTPREFIX"]], i)
    donefile = sprintf("%s.done", outfile)
    if (file.exists(donefile))
        next
    message(sprintf("Fitting cardinality %d", i))
    fit = identifySignatures(burden, i, nmfDecomposition, nrun = B, includeFit = TRUE, seed = seed, .options = sprintf("p%d", cores))
    saveRDS(fit, outfile)
    file.create(donefile)
}

file.create(sprintf("%s.done", env[["OUTPREFIX"]]))
