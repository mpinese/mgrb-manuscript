setwd("C:/Users/mpinese/Downloads/MGRB_revisions/common_var_burden_expanded/")
sb = readRDS("stressboot_dfs.rds")

library(ggplot2)
sb$boot$deriv_cohort = c("aus_from_ukbb" = "UKBB - MGRB", "aus_from_gnomad" = "GnomAD - MGRB")[sb$boot$deriv_cohort]
sb$boot$model = gsub(":[^:]+$", "", sb$boot$model)
ggplot(sb$boot, aes(x = value.rel, y = stat(ncount))) + 
    geom_histogram(binwidth = 2e-5) + facet_grid(deriv_cohort ~ model, scales = "free") + 
    theme_bw() + geom_vline(xintercept = 0, col = "red", lwd = 1.5, lty = "dashed")


ggsave("bootplots_01.svg", ggplot(sb$boot[grepl("^[A-D]", sb$boot$model),], aes(x = value.rel, y = stat(ncount))) + 
    geom_histogram(binwidth = 2e-5) + facet_grid(deriv_cohort ~ model, scales = "free") + 
    theme_bw() + geom_vline(xintercept = 0, col = "red", lwd = 1.5, lty = "dashed") + xlab("PS relative to MGRB"), width = 10, height = 4)

ggsave("bootplots_02.svg", ggplot(sb$boot[grepl("^[E-P]", sb$boot$model),], aes(x = value.rel, y = stat(ncount))) + 
    geom_histogram(binwidth = 2e-5) + facet_grid(deriv_cohort ~ model, scales = "free") + 
    theme_bw() + geom_vline(xintercept = 0, col = "red", lwd = 1.5, lty = "dashed") + xlab("PS relative to MGRB"), width = 10, height = 4)

ggsave("bootplots_03.svg", ggplot(sb$boot[grepl("^[Q-Z]", sb$boot$model),], aes(x = value.rel, y = stat(ncount))) + 
    geom_histogram(binwidth = 2e-5) + facet_grid(deriv_cohort ~ model, scales = "free") + 
    theme_bw() + geom_vline(xintercept = 0, col = "red", lwd = 1.5, lty = "dashed") + xlab("PS relative to MGRB"), width = 6.4, height = 4)


