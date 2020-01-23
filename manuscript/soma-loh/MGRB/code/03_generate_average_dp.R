#!/usr/bin/env Rscript

temp.loci = scan("../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.commonhets.variant_list", sep = ":", what = list(character(0), integer(0), NULL, NULL))
loci = data.frame(chrom = temp.loci[[1]], pos = temp.loci[[2]])
loci = loci[!duplicated(loci),]
loci$key = paste(loci$chrom, loci$pos, sep = ":")

sample_files = list.files("../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.commonhets", "\\.afs\\.bgz$", full.names = TRUE)

sample_files.slices = split(sample_files, 1:28)


makeSliceDepthSumVector = function(slice.sample_files)
{
  all_reldp = matrix(NA, nrow = nrow(loci), ncol = length(slice.sample_files))
  rownames(all_reldp) = loci$key
  colnames(all_reldp) = slice.sample_files

  for (sample_file in slice.sample_files)
  {
    message(sample_file)
    sample_data = scan(sample_file, what = list(character(0), integer(0), integer(0), NULL), quiet = TRUE)
    sample_loci_key = paste(sample_data[[1]], sample_data[[2]], sep = ":")
    sample_dp = sample_data[[3]]
    sample_dp = sample_dp[match(rownames(all_reldp), sample_loci_key)]
    all_reldp[,sample_file] = sample_dp / mean(sample_dp, na.rm = TRUE)
  }

  list(sum = rowSums(all_reldp, na.rm = TRUE), sum2 = rowSums(all_reldp^2, na.rm = TRUE), n = rowSums(!is.na(all_reldp)))
}


library(parallel)

depth_sum_vectors = mclapply(sample_files.slices, makeSliceDepthSumVector, mc.cores = 28)

depth_sum_total = rep(0, nrow(loci))
depth_sum2_total = rep(0, nrow(loci))
depth_count_total = rep(0, nrow(loci))
for (dsv in depth_sum_vectors)
{
  depth_sum_total = depth_sum_total + dsv$sum
  depth_sum2_total = depth_sum2_total + dsv$sum2
  depth_count_total = depth_count_total + dsv$n
}

# mean_reldp and var_reldp are the mean and variance
# of relative depth across loci.
mean_reldp = depth_sum_total/depth_count_total
var_reldp = depth_sum2_total/depth_count_total - mean_reldp^2

# Scale mean_reldp and var_reldp to a mean of one
# across all loci.
scaling_factor = 1 / mean(mean_reldp, na.rm = TRUE)
mean_reldp = mean_reldp * scaling_factor
var_reldp = var_reldp * scaling_factor^2

result = data.frame(chrom = loci$chrom, pos = loci$pos, mean = mean_reldp, var = var_reldp)

saveRDS(result, "../03_dp_stats.rds")
write.table(result, file = "../03_dp_stats.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
