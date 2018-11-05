#!/usr/bin/env Rscript

library(GenomicRanges)
library(ggplot2)
library(plyr)
library(SomaticSignatures)


pdf("../01_highdepth_factorise.pdf", height = 8, width = 8)


# Load data
# Dubbo-1002_1 is an acknowledged bad sample (admixed?)
load("../input/high-depth/06_pass1_variants.rda")
variants = variants[variants$sample != "Dubbo-1002_1",]
# Fix R's sneaky "T" => "TRUE" autoconversion
variants$alt = as.character(variants$alt)
variants$alt[variants$alt == "TRUE"] = "T"


# Passing variants are not failed by either the frequent 
# (BL2) or common SNV (BLC) filters
variants.pass = variants[!(variants$fails_bl2 | variants$fails_blc),]


# Convert the variants.pass data frame to a signatures
# matrix.  Elements of the matrix contain a background-
# corrected count of each signature, per sample.
#
# The matrix is scaled as mutations per megabase.
sigs = ddply(as.data.frame(variants.pass), .(sample, context, alt), nrow)
colnames(sigs)[colnames(sigs) == "V1"] = "count"
temp = ddply(backgrounds, .(sample, context), function(d) sum(d$background))
colnames(temp)[colnames(temp) == "V1"] = "background"
sigs = merge(sigs, temp)
sigs$count.normed = sigs$count / sigs$background * 1e6

temp = expand.grid(base1 = c("A", "C", "G", "T"), base2 = c("C", "T"), base3 = c("A", "C", "G", "T"), alt = c("A", "C", "G", "T"))
temp = temp[as.character(temp$base2) != as.character(temp$alt),]
temp = temp[order(temp$base2, temp$alt, temp$base1, temp$base3),]
all_sigs = sprintf("%s%s%s %s", temp$base1, temp$base2, temp$base3, temp$alt)

sigs_mat = daply(sigs, .(sample, paste(context, alt)), function(d) d$count.normed)
sigs_mat = sigs_mat[,match(all_sigs, colnames(sigs_mat))]
colnames(sigs_mat) = all_sigs
sigs_mat[is.na(sigs_mat)] = 0
image(sigs_mat)
sigs_mat.total = colSums(sigs_mat)
barplot(sigs_mat.total, las=2, family = "mono", cex.axis = 2, col = rep(c("blue", "black", "red", "grey", "green", "pink"), each = 16))

rm(all_sigs)


# Convert the signatures matrix to SomaticSignatures format.
makeSSmotif = function(motif)
{
    sprintf("%s%s %s.%s", substr(motif, 2, 2), substr(motif, 5, 5), substr(motif, 1, 1), substr(motif, 3, 3))
}

motif_mat = t(sigs_mat)
rownames(motif_mat) = makeSSmotif(rownames(motif_mat))
names(dimnames(motif_mat)) = NULL


# Add a small amount of noise for stability
set.seed(1234)
motif_mat = motif_mat + abs(rnorm(prod(dim(motif_mat)), sd = 0.001))

# Factorise
# Programmatic assessment of signature count isn't very reliable here,
# but do it anyway for completeness.  Visual inspection has been more
# fruitful, and indicated that 2 signatures seems to be the best --
# more than that and the canonical deamination signature splits, which
# we know isn't reasonable.
set.seed(1234)
gof_nmf = assessNumberSignatures(motif_mat, 2:15, nReplicates = 5, nrun = 10)
plotNumberSignatures(gof_nmf)
sigs_nmf4 = identifySignatures(motif_mat, 4, nmfDecomposition, nrun = 50)
sigs_nmf3 = identifySignatures(motif_mat, 3, nmfDecomposition, nrun = 50)
sigs_nmf2 = identifySignatures(motif_mat, 2, nmfDecomposition, nrun = 50)
plotSignatureMap(sigs_nmf4)
plotSignatureMap(sigs_nmf3)
plotSignatureMap(sigs_nmf2)
plotSignatures(sigs_nmf4)
plotSignatures(sigs_nmf3)
plotSignatures(sigs_nmf2)
plotSampleMap(sigs_nmf4)
plotSampleMap(sigs_nmf3)
plotSampleMap(sigs_nmf2)
plotSamples(sigs_nmf4)
plotSamples(sigs_nmf3)
plotSamples(sigs_nmf2)


# Generate some plots
# Rescale the NMF weights back to mutations per megabase.
temp = as.data.frame(t(t(sigs_nmf2@samples) * colSums(sigs_nmf2@signatures)))
temp$Sample = gsub("_.*", "", rownames(temp))
temp$Cohort = gsub("-[^_]*", "", rownames(temp))
temp$Burden = colSums(motif_mat)[rownames(temp)]

ggplot(temp, aes(x = Cohort, y = Burden, col = Cohort)) + geom_violin(trim = FALSE) + geom_dotplot(aes(fill = Cohort), binaxis = "y", stackdir = "center", dotsize = 0.25, stackratio = 0.9) + ggtitle("Burden")
ggplot(temp, aes(x = Cohort, y = S1, col = Cohort)) + geom_violin(trim = FALSE) + geom_dotplot(aes(fill = Cohort), binaxis = "y", stackdir = "center", dotsize = 0.25, stackratio = 0.9) + ggtitle("Signature 1")
ggplot(temp, aes(x = Cohort, y = S2, col = Cohort)) + geom_violin(trim = FALSE) + geom_dotplot(aes(fill = Cohort), binaxis = "y", stackdir = "center", dotsize = 0.25, stackratio = 0.9) + ggtitle("Signature 2")

ggplot(temp, aes(x = Cohort, y = Burden, col = Cohort)) + geom_boxplot() + ggtitle("Total mutation burden")
ggplot(temp, aes(x = Cohort, y = S1, col = Cohort)) + geom_boxplot() + ggtitle("Signature 1 burden")
ggplot(temp, aes(x = Cohort, y = S2, col = Cohort)) + geom_boxplot() + ggtitle("Signature 2 burden")


# Tentatively assign signature source to variants.  Derive as:
# Pr(S1|Motif,Sample) = Pr(Motif|S1,Sample)*Pr(S1|Sample) / (Pr(Motif|S1,Sample)*Pr(S1|Sample) + Pr(Motif|S2,Sample)*Pr(S2|Sample))
# And assuming that Pr(Motif|S1,Sample) = Pr(Motif|S1)  (ie that the signatures have identical effect between samples, necessarily)
temp.pr_s1_sample = sigs_nmf2@samples[,1] / rowSums(sigs_nmf2@samples)
temp.pr_s2_sample = 1 - temp.pr_s1_sample
temp.pr_motif_s1 = sigs_nmf2@signatures[,1] / sum(sigs_nmf2@signatures[,1])
temp.pr_motif_s2 = sigs_nmf2@signatures[,2] / sum(sigs_nmf2@signatures[,2])
variants$motif = paste(substr(variants$context, 2, 2), variants$alt, " ", substr(variants$context, 1, 1), ".", substr(variants$context, 3, 3), sep = "")
variants$pr_s1 = temp.pr_motif_s1[variants$motif]*temp.pr_s1_sample[as.vector(variants$sample)] / 
    (temp.pr_motif_s1[variants$motif]*temp.pr_s1_sample[as.vector(variants$sample)] + 
     temp.pr_motif_s2[variants$motif]*temp.pr_s2_sample[as.vector(variants$sample)])



# Export the signatures, loadings, and variants
write.table(sigs_nmf2@samples, "../01_highdepth_loadings.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(sigs_nmf2@signatures, "../01_highdepth_signatures.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(variants, "../01_highdepth_variants.tsv.xz", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)


# Save the workspace
rm(list = grep("^temp", ls(), value = TRUE))
save.image("../01_highdepth_factorise.rda")




# Examine the T>C signature in more detail.  The hypothesis is that
# it's either machine error, or GL variation
#   If in bl2 NOT blc => recurrently mutated loci => machine error
#   If in blc NOT bl2 => sites of common variation => GL variation
temp = expand.grid(base1 = c("A", "C", "G", "T"), base2 = c("C", "T"), base3 = c("A", "C", "G", "T"), alt = c("A", "C", "G", "T"))
temp = temp[as.character(temp$base2) != as.character(temp$alt),]
temp = temp[order(temp$base2, temp$alt, temp$base1, temp$base3),]
all_sigs = sprintf("%s%s%s %s", temp$base1, temp$base2, temp$base3, temp$alt)

sigs.bl2_not_blc = ddply(as.data.frame(variants[variants$fails_bl2 & (!variants$fails_blc),]), .(sample, context, alt), nrow)
sigs.bl2_not_blc = sigs.bl2_not_blc[sigs.bl2_not_blc$sample != "Dubbo-1002_1",]
colnames(sigs.bl2_not_blc)[colnames(sigs.bl2_not_blc) == "V1"] = "count"
temp = ddply(backgrounds, .(sample, context), function(d) sum(d$background))
colnames(temp)[colnames(temp) == "V1"] = "background"
sigs.bl2_not_blc = merge(sigs.bl2_not_blc, temp)
sigs.bl2_not_blc$count.normed = sigs.bl2_not_blc$count / sigs.bl2_not_blc$background * 1e6
sigs_mat.bl2_not_blc = daply(sigs.bl2_not_blc, .(sample, paste(context, alt)), function(d) d$count.normed)
sigs_mat.bl2_not_blc = sigs_mat.bl2_not_blc[,match(all_sigs, colnames(sigs_mat.bl2_not_blc))]
colnames(sigs_mat.bl2_not_blc) = all_sigs
sigs_mat.bl2_not_blc[is.na(sigs_mat.bl2_not_blc)] = 0
sigs_mat.bl2_not_blc = sigs_mat.bl2_not_blc[,order(substr(colnames(sigs_mat.bl2_not_blc), 2, 2), substr(colnames(sigs_mat.bl2_not_blc), 5, 5), substr(colnames(sigs_mat.bl2_not_blc), 1, 1), substr(colnames(sigs_mat.bl2_not_blc), 3, 3))]
image(sigs_mat.bl2_not_blc)
sigs_mat.bl2_not_blc.total = colSums(sigs_mat.bl2_not_blc)

sigs.blc_not_bl2 = ddply(as.data.frame(variants[(!variants$fails_bl2) & variants$fails_blc,]), .(sample, context, alt), nrow)
sigs.blc_not_bl2 = sigs.blc_not_bl2[sigs.blc_not_bl2$sample != "Dubbo-1002_1",]
colnames(sigs.blc_not_bl2)[colnames(sigs.blc_not_bl2) == "V1"] = "count"
temp = ddply(backgrounds, .(sample, context), function(d) sum(d$background))
colnames(temp)[colnames(temp) == "V1"] = "background"
sigs.blc_not_bl2 = merge(sigs.blc_not_bl2, temp)
sigs.blc_not_bl2$count.normed = sigs.blc_not_bl2$count / sigs.blc_not_bl2$background * 1e6
sigs_mat.blc_not_bl2 = daply(sigs.blc_not_bl2, .(sample, paste(context, alt)), function(d) d$count.normed)
sigs_mat.blc_not_bl2 = sigs_mat.blc_not_bl2[,match(all_sigs, colnames(sigs_mat.blc_not_bl2))]
colnames(sigs_mat.blc_not_bl2) = all_sigs
sigs_mat.blc_not_bl2[is.na(sigs_mat.blc_not_bl2)] = 0
sigs_mat.blc_not_bl2 = sigs_mat.blc_not_bl2[,order(substr(colnames(sigs_mat.blc_not_bl2), 2, 2), substr(colnames(sigs_mat.blc_not_bl2), 5, 5), substr(colnames(sigs_mat.blc_not_bl2), 1, 1), substr(colnames(sigs_mat.blc_not_bl2), 3, 3))]
image(sigs_mat.blc_not_bl2)
sigs_mat.blc_not_bl2.total = colSums(sigs_mat.blc_not_bl2)

stopifnot(names(sigs_mat.bl2_not_blc.total) == names(sigs_mat.total))
stopifnot(names(sigs_mat.blc_not_bl2.total) == names(sigs_mat.total))
par(mfrow = c(5, 1))
barplot(sigs_mat.total, las=2, family = "mono", cex.axis = 2, col = rep(c("blue", "black", "red", "grey", "green", "pink"), each = 16), main = "Pass (!BL2 !BLC)")
barplot(sigs_nmf2@signatures[,1], las=2, family = "mono", cex.axis = 2, col = rep(c("blue", "black", "red", "grey", "green", "pink"), each = 16), main = "Pass S1 (!BL2 !BLC S1)")
barplot(sigs_nmf2@signatures[,2], las=2, family = "mono", cex.axis = 2, col = rep(c("blue", "black", "red", "grey", "green", "pink"), each = 16), main = "Pass S2 (!BL2 !BLC S2)")
barplot(sigs_mat.bl2_not_blc.total, las=2, family = "mono", cex.axis = 2, col = rep(c("blue", "black", "red", "grey", "green", "pink"), each = 16), main = "Recurrent (BL2 !BLC)")
barplot(sigs_mat.blc_not_bl2.total, las=2, family = "mono", cex.axis = 2, col = rep(c("blue", "black", "red", "grey", "green", "pink"), each = 16), main = "Common (!BL2 BLC)")
par(mfrow = c(1, 1))

# A bit of an equivocal result -- the odd signature is present in both
# the recurrent and the common sets.  Not sure what that means... unless
# it's a pervasive machine noise signal, and therefore is falsely present
# in the GnomAD database (as common variants), as well as often appears
# in our data (as recurrent 'somatic' variants).


dev.off()

