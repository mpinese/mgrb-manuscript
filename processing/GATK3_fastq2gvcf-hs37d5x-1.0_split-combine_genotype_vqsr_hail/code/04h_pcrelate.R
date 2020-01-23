options(echo = TRUE)
library(SNPRelate)

infile_stem = "../ldpruned/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.ldpruned.plink"
outfile_gds = "../ldpruned/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.ldpruned.gds"
outfile_pcair_rds = "../ldpruned/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.ldpruned.pcair.rds"
outfile_pcrelate_gds = "../ldpruned/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.ldpruned.pcrelate"
outfile_pcrelate_rds = "../ldpruned/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.ldpruned.pcrelate.rds"
outfile_svg = "../ldpruned/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.ldpruned.pcrelate.svg"
outfile_samples_to_drop = "../ldpruned/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.ldpruned.pcrelate.samples_to_drop"

# Convert the PLINK BED to GDS for reading
snpgdsBED2GDS(bed.fn = sprintf("%s.bed", infile_stem), bim.fn = sprintf("%s.bim", infile_stem), fam.fn = sprintf("%s.fam", infile_stem), out.gdsfn = outfile_gds, snpfirstdim = TRUE)

# Convert sample IDs to numeric, as needed by pcrelate below
library(gdsfmt)
genoFile <- openfn.gds(outfile_gds, readonly = FALSE)
sampleid_node = index.gdsn(genoFile, "sample.id")
rename.gdsn(sampleid_node, "sample.name")
samplename_node = sampleid_node
samplename = read.gdsn(samplename_node)
sampleid_node = add.gdsn(genoFile, "sample.id", val = 1:length(samplename))
snpid_node = index.gdsn(genoFile, "snp.id")
genotype_node = index.gdsn(genoFile, "genotype")
moveto.gdsn(sampleid_node, samplename_node, "before")
moveto.gdsn(snpid_node, sampleid_node, "after")
moveto.gdsn(genotype_node, snpid_node, "after")
closefn.gds(genoFile)

# Calculate relatedness estimates using the KING algorithm
genoFile <- snpgdsOpen(outfile_gds)
library(parallel)
KINGmat <- snpgdsIBDKING(genoFile, type = "KING-robust", num.thread = detectCores())
KINGmat2 = KINGmat$kinship
rownames(KINGmat2) = KINGmat$sample.id
colnames(KINGmat2) = KINGmat$sample.id
snpgdsClose(genoFile)

library(GENESIS)
library(GWASTools)

# Read genotypes back in
genoFile <- GdsGenotypeReader(filename = outfile_gds)
genoData <- GenotypeData(genoFile)

# Run PC-AiR to partition the sample
pcair_auto <- pcair(genoData = genoData, kinMat = KINGmat2, divMat = KINGmat2)

saveRDS(list(auto = pcair_auto), outfile_pcair_rds)

summary(pcair_auto)

pcrelate(genoData = genoData, pcMat = pcair_auto$vectors[,1:3], training.set = pcair_auto$unrels, write.to.gds = TRUE, gds.prefix = outfile_pcrelate_gds)

close(genoFile)


# Examine the PCrelate results
library(gdsfmt)
library(GENESIS)

genoFile <- openfn.gds(outfile_gds, readonly = FALSE)
sample_ids = read.gdsn(index.gdsn(genoFile, "sample.id"))
sample_names = read.gdsn(index.gdsn(genoFile, "sample.name"))
closefn.gds(genoFile)

pcrelate <- openfn.gds(sprintf("%s_pcrelate.gds", outfile_pcrelate_gds))
kinship = pcrelateReadKinship(pcrelObj = pcrelate, kin.thresh = 2^(-11/2))
inbreed = pcrelateReadInbreed(pcrelObj = pcrelate, f.thresh = 2^(-11/2))
grm = pcrelateMakeGRM(pcrelObj = pcrelate, scaleKin = 2)
closefn.gds(pcrelate)

kinship$name1 = sample_names[match(kinship$ID1, sample_ids)]
kinship$name2 = sample_names[match(kinship$ID2, sample_ids)]
inbreed$name = sample_names[match(inbreed$ID, sample_ids)]
rownames(grm) = sample_names[match(rownames(grm), sample_ids)]
colnames(grm) = sample_names[match(colnames(grm), sample_ids)]

pcrelate_results = list(kinship = kinship, inbreed = inbreed, grm = grm)

saveRDS(pcrelate_results, outfile_pcrelate_rds)

kinship_MGRB = pcrelate_results$kinship[grepl("^[ABZ]", pcrelate_results$kinship$name1) | grepl("^[ABZ]", pcrelate_results$kinship$name2),]
inbreed_MGRB = pcrelate_results$inbreed[grepl("^[ABZ]", pcrelate_results$inbreed$name),]
grm_MGRB = pcrelate_results$grm[grepl("^[ABZ]", rownames(pcrelate_results$grm)),grepl("^[ABZ]", colnames(pcrelate_results$grm))]


kinship_MGRB$degree = "5+"
for (degree_i in 0:4)
{
    kinship_MGRB$degree[kinship_MGRB$kin >= 2^(-(degree_i+3/2)) & kinship_MGRB$kin < 2^(-(degree_i+1/2))] = as.character(degree_i)
}
kinship_MGRB$degree = ordered(as.character(kinship_MGRB$degree), levels = c(as.character(0:4), "5+"))

library(ggplot2)

temp_segments_1 = data.frame(
    x =      -1, 
    y =      c(0.5,  0.25, 0.125, 0.0625, 0.03125, 0.015625), 
    xend =   c(0,    0.25, 0.5,   0.75,   0.875,   0.9375  ), 
    yend =   c(0.5,  0.25, 0.125, 0.0625, 0.03125, 0.015625),
    degree = c("0",  "1",  "2",   "3",    "4",    "5+"     ))
temp_segments_2 = data.frame(
    x =      c(0,    0.25, 0.5,   0.75,   0.875,   0.9375  ), 
    y =      -1, 
    xend =   c(0,    0.25, 0.5,   0.75,   0.875,   0.9375  ), 
    yend =   c(0.5,  0.25, 0.125, 0.0625, 0.03125, 0.015625),
    degree = c("0", "1",  "2",   "3",    "4",     "5+"     ))
temp_segments = rbind(temp_segments_1, temp_segments_2)

library(ggplot2)
svg(outfile_svg)
ggplot(kinship_MGRB, aes(x = k0, y = kin, col = degree)) + 
    geom_abline(slope = -0.25, intercept = 0.25, linetype = "dashed", colour = "grey") + 
    geom_point() + 
    scale_x_continuous(breaks = c(0, 1/4, 1/2, 3/4, 1)) + scale_y_continuous(breaks = c(1/8, 2/8, 3/8, 4/8)) + coord_fixed(2, xlim = c(0, 1), ylim = c(0, 0.5)) + 
    labs(x = expression(k^(0)), y = "Kinship coefficient", col = "Degree") + theme_bw() + 
    geom_segment(data = temp_segments, aes(x = x, y = y, xend = xend, yend = yend))
dev.off()

head(inbreed_MGRB[order(-inbreed_MGRB$f),], 20)
head(kinship_MGRB[order(-kinship_MGRB$kin),], 50)


# Identify a set of samples to drop to yield an unrelated set.
# Algorithm:
#   If a pair is identical (degree == 0):
#     If the pair is in the same cohort:
#       Drop both
#     Else:
#       Drop one
#   Else if a pair is related (degree <= 2):
#     Drop one
#
# Consider only pairs of two MGRB samples (A* or B* sample IDs) 
# when deciding samples to drop.
#
# When dropping one of two samples, arbitrarily choose the one
# with the lowest ID.

kinship.todrop = kinship_MGRB[grepl("^[AB]", kinship_MGRB$name1) & grepl("^[AB]", kinship_MGRB$name2),]

kinship.todrop.identical = kinship.todrop[kinship.todrop$degree == 0,]
kinship.todrop.identical.within_cohort = kinship.todrop.identical[substr(kinship.todrop.identical$name1, 1, 1) == substr(kinship.todrop.identical$name2, 1, 1),]
kinship.todrop.identical.between_cohort = kinship.todrop.identical[substr(kinship.todrop.identical$name1, 1, 1) != substr(kinship.todrop.identical$name2, 1, 1),]
kinship.todrop.identical.between_cohort$todrop = kinship.todrop.identical.between_cohort$name1
kinship.todrop.identical.between_cohort$todrop[kinship.todrop.identical.between_cohort$ID1 < kinship.todrop.identical.between_cohort$ID2] = kinship.todrop.identical.between_cohort$name2

kinship.todrop.related = kinship.todrop[kinship.todrop$degree > 0 & kinship.todrop$degree <= 2,]
kinship.todrop.related$todrop = kinship.todrop.related$name1
kinship.todrop.related$todrop[kinship.todrop.related$ID1 < kinship.todrop.related$ID2] = kinship.todrop.related$name2

samples_to_drop = c(kinship.todrop.identical.within_cohort$name1, kinship.todrop.identical.within_cohort$name2)
samples_to_drop = c(samples_to_drop, kinship.todrop.identical.between_cohort$todrop)
samples_to_drop = c(samples_to_drop, kinship.todrop.related$todrop)

samples_to_drop = sort(samples_to_drop)

write.table(data.frame(sample = samples_to_drop), file = outfile_samples_to_drop, col.names = FALSE, quote = FALSE, row.names = FALSE)
