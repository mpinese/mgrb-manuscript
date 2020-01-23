#!/usr/bin/env Rscript

stats.external = readRDS("../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.manual_polygenic_scores.GnomAD_NFE_AFs.popstar.externalref.permsummary.rds")
stats.internal = readRDS("../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.manual_polygenic_scores.GnomAD_NFE_AFs.popstar.internalref.permsummary.rds")

stats = rbind(stats.external, stats.internal)
rm(stats.external, stats.internal)
stats$ref_af = factor(c("internal", "external")[stats$external_ref_af + 1])
stats$model = factor(stats$model)
stats = stats[,colnames(stats) != "external_ref_af"]

library(reshape2)
stats = melt(stats, id.vars = c("model", "iter", "seed", "nafbins", "ref_af"), value.name = "statistic")

# Statistic m1 makes no sense for internal reference -- expectation is zero
stats = stats[!(stats$statistic == "m1" & stats$ref_af == "internal"),]

# For now consider only "external" reference data -- their interpretation is more straightforward.
stats = stats[stats$ref_af == "external",]

# We have two Breast cancer scores -- Li and Michailidou.  Michailidou is by far the more
# recent, so exclude Li.
stats = stats[stats$model != "BreastCancer:Li:10.1038/gim.2016.43",]

# Deelen's longevity signature is only based on four SNPs, so any higher-order estimates will
# be problematic indeed.  Exclude these.
stats = stats[!(stats$model == "Longevity:Deelen:10.1093/hmg/ddu139" & stats$variable != "m1"),]

library(plyr)
pvals = ddply(stats, .(model, ref_af, variable), function(d) {
    stat_orig = d$statistic[d$iter == 0]
    stat_perm = d$statistic[d$iter != 0]
    more_extreme = min(sum(stat_orig < stat_perm), sum(stat_orig > stat_perm))
    pval = (more_extreme + 1) / (length(stat_perm) + 1)
    zscore = (stat_orig - mean(stat_perm)) / sd(stat_perm)
    c(p = pval, Z = zscore)
})

pvals$pholm = p.adjust(pvals$p, "holm")

temp = pvals[pvals$pholm < 0.05,]
temp = temp[order(as.integer(gsub("^m", "", temp$variable)), temp$Z),]
temp

write.csv(pvals, file = "../03_depletion_pvals.csv", row.names = FALSE)
write.csv(temp, file = "../03_depletion_pvals_sig.csv", row.names = FALSE)


scores.external = readRDS("../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.manual_polygenic_scores.GnomAD_NFE_AFs.popstar.externalref.permscores.rds")

if (!dir.exists("../03_depletion_plots/"))
    dir.create("../03_depletion_plots/")

for (j in unique(scores.external$model))
{
    Cairo::CairoPNG(sprintf("../03_depletion_plots/%s.external1000.density.png", gsub("[:/]", "_", j)), height = 1500, width = 1500, res = 300)
    temp = scores.external[scores.external$model == j,]
    temp.2 = density(temp$value[temp$iter == 0])
    temp.3 = density(temp$value[temp$iter == 1])
    plot(0 ~ 1, type = "n", xlim = range(temp$value), ylim = c(0, max(max(temp.3$y), max(temp.2$y)))*1.1, xlab = "Polygenic score", ylab = "Density", main = j)
    for (i in 1:max(temp$iter))
    {
        lines(density(temp$value[temp$iter == i]), col = rgb(0, 0, 0, 0.01))
    }
    lines(temp.2, col = rgb(1, 0, 0, 1), lwd = 1.5)
    dev.off()

    Cairo::CairoPNG(sprintf("../03_depletion_plots/%s.external1000.ecdf.png", gsub("[:/]", "_", j)), height = 1500, width = 1500, res = 300)
    temp = scores.external[scores.external$model == j,]
    plot(0 ~ 1, type = "n", xlim = range(temp$value), ylim = c(0, 1), xlab = "Polygenic score", ylab = "Cumulative density", main = j)
    temp.x = seq(min(temp$value), max(temp$value), length.out = 1000)
    for (i in 1:max(temp$iter))
    {
        lines(temp.x, ecdf(temp$value[temp$iter == i])(temp.x), col = rgb(0, 0, 0, 0.01))
    }
    lines(temp.x, ecdf(temp$value[temp$iter == 0])(temp.x), col = rgb(1, 0, 0, 1), lwd = 1.5)
    dev.off()
}

