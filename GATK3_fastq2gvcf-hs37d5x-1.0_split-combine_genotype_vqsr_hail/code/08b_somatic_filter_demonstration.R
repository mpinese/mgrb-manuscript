# Demonstration of the effect of the depth-variable somatic filters.

scenario = expand.grid(depth = 1:100, alpha = c(0.05, 0.1, 0.2))

scenario$ad_lower = qbinom(scenario$alpha/2, scenario$depth, 0.5, lower.tail = TRUE)
scenario$ad_upper = qbinom(scenario$alpha/2, scenario$depth, 0.5, lower.tail = FALSE)

scenario$vaf_lower = scenario$ad_lower / scenario$depth
scenario$vaf_upper = scenario$ad_upper / scenario$depth

library(ggplot2)

# ggplot(scenario, aes(x = depth, ymin = vaf_lower, ymax = vaf_upper, fill = ordered(alpha))) + geom_ribbon() + xlab("Depth") + ylab("VAF")

ggplot(scenario[scenario$alpha == 0.1,], aes(x = depth, ymin = vaf_lower, ymax = vaf_upper)) + geom_ribbon(alpha = 0.5) + xlab("Depth") + ylab("VAF") + theme_bw() + geom_hline(yintercept = c(0.3, 1-0.3)) + geom_vline(xintercept = 10)


