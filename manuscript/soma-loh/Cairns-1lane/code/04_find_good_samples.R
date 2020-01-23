#!/usr/bin/env Rscript

dps = read.table("../03_Cairns-1lane.median_depths.tsv", header = TRUE)

dps$dh = dps$h2 - dps$h1

sel = dps$dh <= 10 & dps$med >= 20 & dps$max <= 100 & dps$med*2+15 >= dps$max

write.table(dps$sample[sel], file = "../04_good_samples.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

