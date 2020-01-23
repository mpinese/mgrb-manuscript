setwd("C:/Users/mpinese/Downloads/MGRB_revisions")

library(openxlsx)
library(svglite)

tbl = read.xlsx("FreeBayes somatic variants in MGRB.xlsx")

tbl.sel = tbl[tbl$InMGRBPhase2 == 1 & tbl$HotspotOrTSGLoss_10pct == 1,]

genes = sort(unique(tbl.sel$SYMBOL))
sids = sort(unique(tbl.sel$sampleID))

gene_samp = 1*sapply(genes, function(gene) sids %in% tbl.sel$sampleID[tbl.sel$SYMBOL == gene])
rownames(gene_samp) = sids

gene_samp = gene_samp[,order(-colSums(gene_samp))]

gene_nsamps = colSums(gene_samp)

xcoords = barplot(gene_nsamps, las=2, ylab = "Number of individuals", ylim = c(0, 50))
text(xcoords, gene_nsamps + 1, label = as.character(gene_nsamps))

svglite("somatic_snv_gene_plot.svg", height = 8, width = 8)
par(mar = c(5, 6, 0, 2)/2)
xcoords = barplot(rev(gene_nsamps[gene_nsamps > 1]), las=1, horiz=TRUE, space = 0.7, col="black", border = NA, xlab = "Number of individuals", xlim = c(0, 50), cex.names=0.9, cex.axis=0.9)
text(rev(gene_nsamps[gene_nsamps > 1]) + 1.5, xcoords, label = rev(as.character(gene_nsamps[gene_nsamps > 1])), cex = 0.9)
dev.off()

gene_nsamps2 = c(gene_nsamps[gene_nsamps > 1], Other = sum(gene_nsamps[gene_nsamps <= 1]))
xcoords = barplot(gene_nsamps2, las=2, ylab = "Number of individuals", ylim = c(0, 50))
text(xcoords, gene_nsamps2 + 1, label = as.character(gene_nsamps2))

