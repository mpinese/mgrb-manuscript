#!/usr/bin/env Rscript

library(data.table)

saveRDS(fread("../MGRB_phase2.dupmarked.realigned.recalibrated.cgc.freebayes.split.vtnorm.merged.tsv"), "../MGRB_phase2.dupmarked.realigned.recalibrated.cgc.freebayes.split.vtnorm.merged.rds")
