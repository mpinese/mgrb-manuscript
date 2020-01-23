library(ggplot2)
library(plyr)
library(SomaticSignatures)

bestfits = readRDS("../05_combined_factorise_bestfits.rds")
factsum = readRDS("../05_combined_factorise_summary.rds")

factsum.best = ddply(factsum, .(k, mode), function(d) d[which.max(d$evar),])

# Use the age grouped factorization (with lower noise) to establish
# a good cardinality.

# Use the age grouped factorization (with lower noise) to establish
# a good cardinality.
ggplot(factsum.best[factsum.best$mode == "agegrouped",], aes(x = as.integer(k), y = evar)) + geom_point() + geom_line() +
    theme_bw() + xlab("NMF factorisation cardinality") + ylab("Explained variance") + scale_x_continuous(breaks = 0:12, minor_breaks = NULL) +
    coord_cartesian(xlim = c(0, 12), ylim = c(0.8, 1))
# It's a bit of a tricky call but considering the (hidden) points
# at k = 0, 1, k=3 seems like the point at which the evar improvement
# becomes linear.

# Plot signatures for the single fit at k=3 as this is what's used
# later in the paper.
plotSignatures(bestfits[["single.3"]])
# For reference also look at the agegrouped fit for k=3 to establish
# that the signatures are roughly equivalent.
plotSignatures(bestfits[["agegrouped.3"]])

# A bit less noise in agegrouped, as expected, but definitely concordant.


