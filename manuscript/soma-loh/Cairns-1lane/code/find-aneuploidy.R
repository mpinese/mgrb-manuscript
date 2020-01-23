#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(docopt))
suppressPackageStartupMessages(library(tvd))
suppressPackageStartupMessages(library(plyr))


loadTrace = function(path, vaf_lt = 0.05, vaf_ht = 0.95, lambda_vaf = 2, lambda_dp = 20, dpcalib = NULL, plot_diag = FALSE, ...)
{
    if (is.null(dpcalib))
    {
        dpcalib = data.frame(chrom = character(0), pos = integer(0), dpcalib = numeric(0))
    }
    temp = read.table(path, sep = "\t", stringsAsFactors = FALSE, header = FALSE, col.names = c("chrom", "pos", "dp", "ad"))
    temp = merge(temp, dpcalib, all.x = TRUE, all.y = FALSE)
    temp$dpcalib[is.na(temp$dpcalib)] = 1.0

    temp = temp[order(as.integer(temp$chrom), temp$pos),]

    temp$vaf = temp$ad / temp$dp
    if (plot_diag)
    {
        hist(temp$vaf, xlim = c(0, 1), breaks = 50, xlab = "VAF", ...)
        abline(v = c(vaf_lt, vaf_ht), col = "red", lwd = 2)
    }
    temp = temp[temp$vaf >= vaf_lt & temp$vaf <= vaf_ht,]
    temp$fvaf = pmin(temp$vaf, 1 - temp$vaf)
    temp$sfvaf = tvd1d(temp$fvaf, lambda = lambda_vaf)
    temp$cdp = temp$dp / temp$dpcalib
    temp$rcdp = temp$cdp / median(temp$cdp)
    temp$srcdp = tvd1d(temp$rcdp, lambda = lambda_dp)
    temp
}


findAneuploidy = function(trace, mad_factor = 10, min_conseq = 5)
{
    threshold = median(trace$sfvaf) - 10*mad(trace$sfvaf)
    trace$aneu = trace$sfvaf < threshold
    
    # Find runs in aneuploidy calls.  Do per-chromosome to guarantee no
    # inter-chromosome overruns.
    runs = ddply(trace, .(chrom), function(d) {
        rl = rle(d$aneu)
        rlc = cumsum(rl$lengths)

        # State machine because I'm lazy
        # State value       Description
        # 0                 Init
        # 1                 In non-aneuploid section
        # 2                 In aneuploid section
        state = 0
        run_starts = c()
        run_ends = c()
        for (i in 1:length(rl$lengths))
        {
            v = rl$values[i]
            w = rl$lengths[i]
            s = rlc[i]

            if (state == 0) 
            {
                if (v == TRUE)  {   state = 2; run_starts = c(1)                        }
                else            {   state = 1                                           }
            }
            else if (state == 1)
            {
                if (v == TRUE)  {   state = 2; run_starts = c(run_starts, s - w + 1)    }
            }
            else if (state == 2)
            {
                if (v == FALSE) {   state = 1; run_ends = c(run_ends, s - w)            }
            }
        }
        if (state == 2) {
            run_ends = c(run_ends, nrow(d))
        }

        if (length(run_starts) == 0)
        {
            return (NULL)
        }

        # Runs will span d[run_start:run_end,], inclusively
        runs = data.frame(start_index = run_starts, end_index = run_ends)
        runs$chrom = d$chrom[[1]]
        runs$start_pos = d$pos[runs$start_index]
        runs$end_pos = d$pos[runs$end_index]
        runs$n_loci = runs$end_index - runs$start_index + 1
        runs$min_size = runs$end_pos - runs$start_pos + 1
        runs$median_sfvaf = mapply(function(start, end) median(d$sfvaf[start:end]), run_starts, run_ends)
        runs$median_rcdp = mapply(function(start, end) median(d$rcdp[start:end]), run_starts, run_ends)
        runs$threshold = threshold

        # Emit only these runs.
        runs
    })

    runs = runs[runs$n_loci >= min_conseq,]
    if (nrow(runs) == 0)
    {
        runs = data.frame(start_index = integer(0), end_index = integer(0), chrom = character(0), start_pos = integer(0), end_pos = integer(0),
            n_loci = integer(0), min_size = integer(0), median_sfvaf = numeric(0), median_rcdp = numeric(0), threshold = numeric(0))
    }

    runs
}


plotRaw = function(trace, jitter_vaf = TRUE, ...)
{
    npoints = sum(trace$vaf >= 0.2 & trace$vaf < 0.8)
    alpha = max(1/254, 1e5/npoints)
    if (jitter_vaf == TRUE)
    {
        jitter_amnt = 0.25/sqrt(mean(trace$dp))
        vaf = jitter(trace$vaf, amount = jitter_amnt)
    }
    else
    {
        vaf = trace$vaf
    }

    plot(vaf, pch = ".", col = rgb(0, 0, 0, alpha), ylim = c(0, 1), xlab = "Position", ylab = "VAF", ...)
    # Inspection suggests that vaf_jittered is very close to normal in the null.
    # qqnorm(sample(vaf, 1000))
}


plotTrace = function(trace, ...)
{
    plot(2*(0.5-trace$sfvaf), pch = ".", ylim = c(0, 1.2), type = "l", lwd = 2, xlab = "Position", ylab = "LoH Score", ...)
    abline(h = 2*(0.5-(median(trace$sfvaf) - 10*mad(trace$sfvaf))), col = "red")
    lines(trace$srcdp, pch = ".", type = "l", lwd = 2, col = "blue")
    abline(v = which(trace$chrom[-1] != trace$chrom[-length(trace$chrom)]), col = rgb(0, 0, 0, 0.2), lwd = 2)
    abline(h = 1, col = "blue", lty = "dotted")
}



VERSION = "v0.0.2  13 Nov 2017"

sprintf('Find regions of subclonal aneuploidy in NGS count data.

Usage:
  find-aneuploidy.R [--lvt=<LVT>] [--uvt=<UVT>] [--lambdav=<LV>] [--lambdad=<LD>] [--madk=<K>] [--minrun=<N>] [--dpcalib=<DPC>] [--diag=<PDF>] <infile> <sampleid> <outfile>
  find-aneuploidy.R -h | --help
  find-aneuploidy.R --version

Options:
 -h --help        Show this message.
 --version        Show version.
 --lvt=<LVT>      Lower VAF threshold [default: 0.07].
 --uvt=<UVT>      Upper VAF threshold [default: 0.93].
 --lambdav=<LV>   Lambda smoothing factor for VAF trace [default: 5.0].
 --lambdad=<LD>   Lambda smoothing factor for depth trace [default: 25.0].
 --madk=<K>       VAF MAD scaling factor for aneuploidy threshold [default: 10.0].
 --minrun=<N>     Minimum number of consecutive aneuploid loci for a call [default: 5].
 --dpcalib=<DPC>  Path to input depth calibration file, tsv format with header, columns chrom, pos, dpcal.
 --diag=<PDF>     Path to PDF of diagnostic plots (if not specified, plots are not generated).

%s
Mark Pinese  <m.pinese@garvan.org.au>', VERSION) -> doc

opts = docopt(doc)

if (opts$version)
{
    cat(VERSION)
    cat("\n")
    quit(save = "no")
}

diag_mode = FALSE
if (!is.null(opts$diag))
{
    diag_mode = TRUE
    pdf(opts$diag, height = 8, width = 8)
}

vaf_lt = as.numeric(opts$lvt)
vaf_ht = as.numeric(opts$uvt)
lambda_vaf = as.numeric(opts$lambdav)
lambda_dp = as.numeric(opts$lambdad)
mad_factor = as.numeric(opts$madk)
min_conseq = as.integer(opts$minrun)

stopifnot(vaf_lt >= 0 && vaf_lt <= 1)
stopifnot(vaf_ht >= 0 && vaf_ht <= 1)
stopifnot(vaf_lt < vaf_ht)
stopifnot(lambda_vaf > 0)
stopifnot(lambda_dp > 0)
stopifnot(mad_factor > 0)
stopifnot(min_conseq >= 1)

dpcalib = data.frame(chrom = character(0), pos = integer(0), dpcalib = numeric(0))
if (!is.null(opts$dpcalib))
{
    dpcalib = read.table(opts$dpcalib, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
    stopifnot(colnames(dpcalib[1:2]) == c("chrom", "pos"))
    stopifnot(ncol(dpcalib) == 3)
    colnames(dpcalib)[3] = "dpcalib"
}

trace = loadTrace(opts$infile, 
    vaf_lt = vaf_lt, vaf_ht = vaf_ht, 
    lambda_vaf = lambda_vaf, lambda_dp = lambda_dp,
    dpcalib = dpcalib,
    plot_diag = diag_mode, main = sprintf("%s raw VAF", opts$sampleid))

if (diag_mode)
{
    hist(trace$vaf, breaks = 100, main = sprintf("%s filtered VAF", opts$sampleid), xlab = "VAF")
    plotRaw(trace, main = sprintf("%s raw trace", opts$sampleid), jitter_vaf = TRUE)
    plotTrace(trace, main = sprintf("%s smoothed trace", opts$sampleid))
}

aneu = findAneuploidy(trace, mad_factor = mad_factor, min_conseq = min_conseq)

if (nrow(aneu) > 0)
{
    aneu = cbind(sample = opts$sampleid, infile = opts$infile, vaf_lt = vaf_lt, vaf_ht = vaf_ht, lambda_vaf = lambda_vaf, lambda_dp = lambda_dp, mad_factor = mad_factor, min_conseq = min_conseq, version = VERSION, aneu)
} else {
    aneu = cbind(sample = character(0), infile = character(0), vaf_lt = numeric(0), vaf_ht = numeric(0), lambda_vaf = numeric(0), lambda_dp = numeric(0), mad_factor = numeric(0), min_conseq = integer(0), version = character(0), aneu)
}

write.table(aneu, file = opts$outfile, row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")



