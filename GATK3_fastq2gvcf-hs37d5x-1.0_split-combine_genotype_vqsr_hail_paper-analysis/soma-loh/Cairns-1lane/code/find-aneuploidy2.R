#!/usr/bin/env Rscript

library(plyr)
library(mgcv)
library(tvd)


# Perform per-sample depth anomaly fitting to realise the model:
#    lambda_ij = D_i*(k_ij*a_j*ac_ij)
# => log(lambda_ij) = log(D_i) + log(k_ij / 2) + log(a_j) + log(ac_ij)
# Where
#   lambda_ij is the expected depth in sample i at locus j
#   D_i       is the global depth of sample i
#   k_ij      is the copy number (diploid = 2) of sample i at locus j (might be fractional due to subclonality)
#   a_j       is the affinity anomaly at locus j (equal to affinity$log2_reldp)
#   ac_ij     is the depth anomaly at locus j for sample i (to be fit in this code)
#
# ac_ij will combine sample-specific GC signals as well as sample-specifc
# affinity signal not captured by the first-order global correction in a_j.
#
# So the goal is to fit ac_ij given d_ij ~ Pois(lambda_ij), D_i, and a_j,
# and assuming k_ij = 2.
#
# From the above, assuming k_ij = 2:
#     log(lambda_ij) = [log(D_i) + log(a_j)] + log(ac_ij)
#  => log(lambda_ij) - [log(D_i) + log(a_j)] ~ s(a_j) + s(gc_j)
# Where
#   gc_j      is the gc content at locus j.
#
# We further split the s(gc_j) into multiple terms based on differing
# window sizes for GC content averaging.  PCA is used to remove 
# co-linearity.
#
# The final output from the fit is a sample-specific estimate of
# log(lambda_ij) - [log(D_i) + log(a_j)] for each locus j.  To this
# estimate we add log(D_i) and log(a_j) to yield the final estimate
# of log(lambda_ij) if k_ij = 2.  This is the output stored in 
# log_lambda_hat
#
# Perform a basic fit for now, setting lambda_ij = d_ij.  The
# large data size should make up for this liberty with the error model.
fitDepthAnomaly = function(sample_data, affinity, gc)
{
    data = merge(sample_data, affinity, by = c("chrom", "pos"), all = FALSE)
    data = merge(data, gc, by = c("chrom", "pos"), all = FALSE)

    logD = log2(mean(data$dp, trim = 0.1))
    logA = log2(data$reldp)

    #gc_pca = prcomp(as.matrix(data[,grepl("^gc[0-9]+$", colnames(data))]), center = TRUE, scale = TRUE)
    gc_pca = prcomp(as.matrix(data[,c("gc100", "gc200", "gc400", "gc600", "gc800")]), center = TRUE, scale = TRUE)

    anomaly = log2(data$dp) - logD - logA
    fit = gam(anomaly ~ s(logA) + s(gc_pca$x[,1]) + s(gc_pca$x[,2]) + s(gc_pca$x[,3]) + s(gc_pca$x[,4]) + s(gc_pca$x[,5]))
    print(anova(fit))
    anomaly_pred = predict(fit)

    data$log_lambda_hat = logD + logA + anomaly_pred

    data
}


testWindowDifference = function(data1, data2)
{
    ad1 = data1$ad
    ad2 = data2$ad
    rd1 = data1$dp - ad1
    rd2 = data2$dp - ad2
    anomaly1 = log2(c(rd1, ad1)) - rep(data1$log_lambda_hat, 2)
    anomaly2 = log2(c(rd2, ad2)) - rep(data2$log_lambda_hat, 2)
    suppressWarnings(ks.test(anomaly1, anomaly2)$p.value)
}


# First split xsomes into small windows.
# Progressively merge by KS test.
# Calculate fits based on windowed regions -- direct optim.
findWindows = function(data, initial_size = 100, alpha = 0.01)
{
    # Define windows
    data = data[order(data$chrom, data$pos),]
    data$original_index = 1:nrow(data)
    windows = ddply(data, .(chrom), function(cd) {
        cd$window_id = paste(cd$chrom, floor(cd$original_index/initial_size)+1, sep = ":")
        windows = data.frame(window_id = unique(cd$window_id))
        windows$start_index = tapply(cd$original_index, cd$window_id, min)[windows$window_id]
        windows$end_index = tapply(cd$original_index, cd$window_id, max)[windows$window_id]
        windows
    })

    # Calculate initial split/merge P-values
    windows$pr_vs_prev = c(-1, mapply(function(start1, end1, start2, end2) testWindowDifference(data[start1:end1,], data[start2:end2,]),
        windows$start[-nrow(windows)], windows$end[-nrow(windows)], windows$start[-1], windows$end[-1]))
    # Chromosomes are never merged
    windows$pr_vs_prev[(c(windows$chrom, "") != c("", windows$chrom))[1:nrow(windows)]] = -1

    # Progressive merging of windows
    ntests = sum(windows$pr_vs_prev >= 0)
    while (TRUE)
    {
        # Find the highest P-value
        pmax_index = which.max(windows$pr_vs_prev)
        pmax = windows$pr_vs_prev[pmax_index]

        # Check for termination
        if (pmax * ntests < alpha)
            break

        # Merge windows pmax_index and pmax_index-1
        # message(sprintf("Merging indices %d (%s:%d-%d) and %d (%s:%d-%d), P-value %f", 
        #     pmax_index-1, windows$chrom[pmax_index - 1], windows$start_index[pmax_index - 1], windows$end_index[pmax_index - 1], 
        #     pmax_index, windows$chrom[pmax_index - 1], windows$start_index[pmax_index], windows$end_index[pmax_index], 
        #     pmax))
        windows$end_index[pmax_index - 1] = windows$end_index[pmax_index]
        windows = windows[-pmax_index,]

        # Recalculate the split/merge P-values for either side of the new window
        if (pmax_index - 2 > 0)
        {
            # Upstream is a valid window.
            if (windows$chrom[pmax_index - 2] != windows$chrom[pmax_index - 1])
                windows$pr_vs_prev[pmax_index - 1] = -1    # Upstream is a different chrom
            else
                windows$pr_vs_prev[pmax_index - 1] = testWindowDifference(data[windows$start[pmax_index-2]:windows$end[pmax_index-2],], data[windows$start[pmax_index-1]:windows$end[pmax_index-1],])
            ntests = ntests + 1
        }
        if (pmax_index - 1 != nrow(windows))
        {
            # Downstream is a valid window
            if (windows$chrom[pmax_index] != windows$chrom[pmax_index - 1])
                windows$pr_vs_prev[pmax_index] = -1        # Downstream is a different chrom
            else
                windows$pr_vs_prev[pmax_index] = testWindowDifference(data[windows$start[pmax_index-1]:windows$end[pmax_index-1],], data[windows$start[pmax_index]:windows$end[pmax_index],])
            ntests = ntests + 1
        }
    }

    # Add coordinates and return
    windows$start_pos = data$pos[windows$start_index]
    windows$end_pos = data$pos[windows$end_index]
    return(windows)
}


llik = function(dp1, dp2, f, k1, k2, lambdahat, w = NULL)
{
    if (is.null(w))
        w = rep(1, length(lambdahat))

    lambda1 = (k1*f + (1-f))/2 * lambdahat
    lambda2 = (k2*f + (1-f))/2 * lambdahat
    d12 = dpois(dp1, lambda1)*dpois(dp2, lambda2)
    d21 = dpois(dp1, lambda2)*dpois(dp2, lambda1)
    d = 0.5*d12 + 0.5*d21
    sum(log(d)*w)
}


fitWindow.kf = function(window_data, ks, f, w = NULL)
{
    # Find the most likely combination of k in the ks, given f,
    # for the data in window_data.  ks is a 2-column matrix with
    # integer columns k1, k2.  Return as list(k1, k2, llik)
    dp1 = window_data$dp - window_data$ad
    dp2 = window_data$ad
    lambdahat = 2^window_data$log_lambda_hat
    lliks = apply(ks, 1, function(k) llik(dp1, dp2, f, k[1], k[2], lambdahat, w))
    besti = which.max(lliks)
    list(k1 = ks[besti,1], k2 = ks[besti,2], f = f, llik = lliks[besti])
}


fitWindows.maxkf = function(data, windows, maxk, f, w = NULL, mink = 0)
{
    # For each window find the most likely combination of 
    # k1, k2 in mink..maxk, given f.  Returned a list of fits,
    # being a data frame of window fits, and npar, being
    # the estimated number of parameters required to encode
    # the k state among the mink..maxk, mink..maxk possibilities.
    # Note that maxk >= 1.
    stopifnot(maxk >= 1 && mink <= maxk)
    ks = expand.grid(k1 = mink:maxk, k2 = mink:maxk)
    ks = ks[ks$k1 >= ks$k2,]
    window_fits = ddply(windows, colnames(windows), function(d) cbind(d, fit = fitWindow.kf(data[d$start_index:d$end_index,], ks, f, w)))
    window_fits
}


findLocalMaxima = function(x)
{
    which(diff(diff(c(-Inf, x, -Inf)) < 0) == 1)
}


fitWindows.maxk = function(data, windows, maxk, w = NULL, plot = FALSE)
{
    # Given maxk, search over f in (0, 1) to find the most
    # likely value.
    objective = function(f) sum(fitWindows.maxkf(data, windows, maxk, f, w)$fit.llik)

    # Assess f at fairly dense equally-spaced points.  Place a slight bias
    # towards low f.
    grid_pow = 0.05
    test.f = rev(seq(0.99^(-grid_pow), 0.05^(-grid_pow), length.out = maxk*10))^(-1/grid_pow)
    test.ll = laply(test.f, objective)
    local_maxima.i = findLocalMaxima(test.ll)

    # Polish each local maximum
    polished_maxima = sapply(local_maxima.i, function(i) {
        lower = test.f[max(1, i - 1)]
        upper = test.f[min(length(test.f), i + 1)]
        optimise(objective, lower = lower, upper = upper, maximum = TRUE)
    })
    global_optimum = polished_maxima[,which.max(polished_maxima["objective",])]

    if (plot)
    {
        plot(test.ll ~ test.f)
        points(polished_maxima["maximum",], polished_maxima["objective",], col = "red", pch = 19)
        points(global_optimum$maximum, global_optimum$objective, col = "red", pch = 1, cex = 2)
    }

    list(
        optimum = data.frame(f = global_optimum$maximum, maxk = maxk, ll = global_optimum$objective), 
        search = data.frame(f = c(test.f, unlist(polished_maxima["maximum",])), maxk = maxk, ll = c(test.ll, unlist(polished_maxima["objective",])))
    )
}


mergeWindowFits = function(window_fits)
{
    # Merge consecutive window fits if they are on the same chromosome
    # and have matching k1, k2, and f values.  Return the merged
    # fits as a data frame.
    if (nrow(window_fits) <= 1)
        return(window_fits)
    window_fits = window_fits[order(window_fits$chrom, window_fits$start_index),]

    merged = window_fits[FALSE,]
    merged_window.key = c(window_fits$chrom[1], window_fits$fit.k1[1], window_fits$fit.k2[1], window_fits$fit.f[1])
    merged_window.start_i = 1

    for (i in 2:nrow(window_fits))
    {
        current_window.key = c(window_fits$chrom[i], window_fits$fit.k1[i], window_fits$fit.k2[i], window_fits$fit.f[i])

        if (any(current_window.key != merged_window.key))
        {
            # A new window.  Add the old one to the merged data.
            merged_window.end_i = i - 1
            merged_window = window_fits[merged_window.start_i,]
            merged_window$end_index = window_fits$end_index[merged_window.end_i]
            merged_window$fit.llik = sum(window_fits$fit.llik[merged_window.start_i:merged_window.end_i])
            merged = rbind(merged, merged_window)

            # Start the new window
            merged_window.start_i = i
            merged_window.key = current_window.key
        }
    }

    # Add the last merged window
    merged_window.end_i = i
    merged_window = window_fits[merged_window.start_i,]
    merged_window$end_index = window_fits$end_index[merged_window.end_i]
    merged_window$fit.llik = sum(window_fits$fit.llik[merged_window.start_i:merged_window.end_i])
    merged = rbind(merged, merged_window)

    # Return without pr_vs_prev, which is invalidated by the merging
    merged[,colnames(merged) != "pr_vs_prev"]
}


fitWindows = function(data, max_maxk, w = NULL, window.initial_size = 100, window.alpha = 0.01)
{
    windows = findWindows(data, initial_size = window.initial_size, alpha = window.alpha)
    model_search = llply(1:max_maxk, function(maxk) fitWindows.maxk(data, windows, maxk, w))
    models = ldply(model_search, function(model_search_result) model_search_result$optimum)
    models$p = ceiling(log2((models$maxk+1)*(models$maxk+2)/2) + 1)
    null_loglik = sum(ddply(windows, colnames(windows), function(d) fitWindow.kf(data[d$start_index:d$end_index,], cbind(k1 = 1, k2 = 1), f = 0, w)$llik)$V1)
    models = rbind(models, c(f = 0, maxk = NA, ll = null_loglik, p = 0))
    models$BIC = -2*models$ll + log(nrow(data))*models$p

    besti = which.min(models$BIC)

    if (is.na(models$maxk[besti]))
        fit = fitWindows.maxkf(data, windows, maxk = 1, f = 0, w = w, mink = 1)
    else
        fit = fitWindows.maxkf(data, windows, maxk = models$maxk[besti], f = models$f[besti], w = w, mink = 0)
    fit = fit[order(fit$chrom, fit$start_index),]

    list(models = models, model_search = ldply(model_search, function(model_search_result) model_search_result$search), fit.orig = fit, fit = mergeWindowFits(fit))
}


plotWindowFit = function(data, fit)
{
    fit$kf1 = fit$fit.f*fit$fit.k1/2 + (1-fit$fit.f)/2
    fit$kf2 = fit$fit.f*fit$fit.k2/2 + (1-fit$fit.f)/2
    fit$d = log2(fit$kf1 + fit$kf2)
    fit$v1 = fit$kf1 / (fit$kf1 + fit$kf2)
    fit$v2 = fit$kf2 / (fit$kf1 + fit$kf2)
    fit$colus = fit$fit.llik / (fit$end_index - fit$start_index + 1)
    fit$col = pmax(0, pmin(1, (fit$colus - (-10)) / 5))
    fit$coldisc = as.integer(ceiling(fit$col * 100))
    fit$coldisc[fit$coldisc == 0] = 1
    pal = rainbow(100, end = 0.33)

    chrom_boundaries = which(data$chrom[-1] != data$chrom[-nrow(data)])
    chrom_midpoints = (c(0, chrom_boundaries) + c(chrom_boundaries, nrow(data))) / 2
    names(chrom_midpoints) = data$chrom[!duplicated(data$chrom)]

    par(mfrow = c(2, 1))
    plot(log2(data$dp) - data$log_lambda_hat, pch = ".", col = rgb(0, 0, 0, 0.25), xlab = "Genomic position", ylab = "Depth anomaly (log2)", main = "Depth", xaxt = "n")
    axis(1, at = chrom_midpoints, labels = names(chrom_midpoints))
#    lines(runmed(log2(data$dp) - data$log_lambda_hat, k = 151), col = "red", lwd = 1)
    lines(tvd1d(log2(data$dp) - data$log_lambda_hat, lambda = 5), col = "red", lwd = 2)
    abline(v = chrom_boundaries, col = "blue", lwd = 2)
    abline(v = fit$start_index, col = "red", lty = "dotted")
    segments(fit$start_index, fit$d, fit$end_index, fit$d, col = pal[fit$coldisc], lwd = 3)

    plot(data$ad / jitter(data$dp, amount = 0.5), pch = ".", ylim = c(0, 1), col = rgb(0, 0, 0, 0.25), xlab = "Genomic position", ylab = "VAF", main = "Allele frequency", xaxt = "n")
    axis(1, at = chrom_midpoints, labels = names(chrom_midpoints))
    abline(v = chrom_boundaries, col = "blue", lwd = 2)
    abline(v = fit$start_index, col = "red", lty = "dotted")
    segments(fit$start_index, fit$v1, fit$end_index, fit$v1, col = pal[fit$coldisc], lwd = 3)
    segments(fit$start_index, fit$v2, fit$end_index, fit$v2, col = pal[fit$coldisc], lwd = 3)
    par(mfrow = c(1, 1))
}


plotRawData = function(data)
{
    chrom_boundaries = which(data$chrom[-1] != data$chrom[-nrow(data)])
    chrom_midpoints = (c(0, chrom_boundaries) + c(chrom_boundaries, nrow(data))) / 2
    names(chrom_midpoints) = data$chrom[!duplicated(data$chrom)]

    par(mfrow = c(2, 1))
    plot(jitter(data$dp, amount = 0.5), pch = ".", col = rgb(0, 0, 0, 0.25), xlab = "Genomic position", ylab = "Depth", main = "Depth", xaxt = "n")
    axis(1, at = chrom_midpoints, labels = names(chrom_midpoints))
    abline(v = chrom_boundaries, col = "blue", lwd = 2)

    plot(data$ad / jitter(data$dp, amount = 0.5), pch = ".", ylim = c(0, 1), col = rgb(0, 0, 0, 0.25), xlab = "Genomic position", ylab = "VAF", main = "Allele frequency", xaxt = "n")
    axis(1, at = chrom_midpoints, labels = names(chrom_midpoints))
    abline(v = chrom_boundaries, col = "blue", lwd = 2)
    par(mfrow = c(1, 1))
}





args = commandArgs(trailingOnly = TRUE)

infile.dp = args[1]
infile.gc = args[2]
outfile.plots = args[3]
infile.afs = args[4]
outfile.aneu = args[5]

message("Loading affinities")
affinity = read.table(infile.dp, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
message("Loading GC")
gc = read.table(infile.gc, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Drop loci with extreme GC
gc = gc[gc$gc100 >= 0.3 & gc$gc100 <= 0.6,]

# Drop loci with extreme affinity
temp.affinity_bounds = quantile(affinity$reldp, c(0.1, 0.9))
affinity = affinity[affinity$reldp >= temp.affinity_bounds[1] & affinity$reldp < temp.affinity_bounds[2],]

# Load afs
message("Loading AFs")
data = read.table(infile.afs, header = FALSE, sep = "\t", stringsAsFactors = FALSE, col.names = c("chrom", "pos", "dp", "ad"))

pdf(outfile.plots, width = 16, height = 8)

# Drop poorly-covered and nearly fixed loci
# Note this filter will need to be adjusted to suit the data source
message("Filtering AFs")
par(mfrow = c(1, 2))
hist(data$ad / data$dp, xlab = "VAF", main = "Pre-filter", breaks = 100, xlim = c(0, 1))
abline(v = c(0.1, 0.9), col = "red")
temp.filter.dp = data$dp >= 5
temp.filter.vaf = data$ad / data$dp >= 0.1 & data$ad / data$dp <= 0.9
message(sprintf("%.2f%% of loci fail depth filter", mean(!temp.filter.dp)*100))
message(sprintf("%.2f%% of loci fail vaf filter", mean(!temp.filter.vaf)*100))
message(sprintf("%.2f%% of loci fail at least one filter", mean(!(temp.filter.dp & temp.filter.vaf))*100))
data = data[temp.filter.dp & temp.filter.vaf,]
hist(data$ad / data$dp, xlab = "VAF", main = "Post-filter", breaks = 100, xlim = c(0, 1))
abline(v = c(0.1, 0.9), col = "red")
par(mfrow = c(1, 1))

# GC and affinity correction
message("Fitting depth anomaly")
data = fitDepthAnomaly(data, affinity, gc)

# Ensure correct ordering
data = data[order(data$chrom, data$pos),]

# Fit
message("Fitting aneuploidy")
fit = fitWindows(data, 4)

# Diagnostic plots
message("Generating plots")
plotRawData(data)
plotWindowFit(data, fit$fit)
dev.off()

# Save results
message("Writing results")
saveRDS(list(data = data, fit = fit, params = list(infile.dp = infile.dp, infile.gc = infile.gc, outfile.plots = outfile.plots, infile.afs = infile.afs, outfile.aneu = outfile.aneu)), outfile.aneu)

message("Done.")
