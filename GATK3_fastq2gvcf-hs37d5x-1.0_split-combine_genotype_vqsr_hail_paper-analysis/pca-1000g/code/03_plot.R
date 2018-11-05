#!/usr/bin/env Rscript

library(ggplot2)

scores = read.table("../02_MGRB_1000G_pca.scores.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

dir.create("../03_plots/", showWarnings = FALSE)

for (dest in c("svg", "png"))
{
    message(dest)
    if (dest == "svg")
        svg("../03_plots/scores_%02d.svg", height = 8, width = 8)
    else if (dest == "png")
        png("../03_plots/scores_%02d.png", height = 800, width = 800)
    print(ggplot(scores, aes(x = PC1, y = PC2, colour = superPopulation)) + geom_point() + theme_bw())
    print(ggplot(scores, aes(x = PC1, y = PC3, colour = superPopulation)) + geom_point() + theme_bw())
    print(ggplot(scores, aes(x = PC1, y = PC4, colour = superPopulation)) + geom_point() + theme_bw())
    print(ggplot(scores, aes(x = PC2, y = PC3, colour = superPopulation)) + geom_point() + theme_bw())
    print(ggplot(scores, aes(x = PC2, y = PC4, colour = superPopulation)) + geom_point() + theme_bw())
    print(ggplot(scores, aes(x = PC3, y = PC4, colour = superPopulation)) + geom_point() + theme_bw())
    dev.off()
}
