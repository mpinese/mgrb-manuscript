options(echo = TRUE)

# Load and merge MGRB and GnomAD data
infile.mgrb = "../02_MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.all.freqs.minimal.vep.filtered.tsv"
infile.gnomad = "../02_gnomad.genomes.r2.0.1.sites.autosomes.split.all.freqs.minimal.vep.filtered.tsv"

vars.mgrb = read.table(infile.mgrb, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
vars.gnomad = read.table(infile.gnomad, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)

temp = strsplit(vars.mgrb$Uploaded_variation, "_")
vars.mgrb$variant = sapply(temp, function(s) s[1])
vars.mgrb$mgrb_freqs = sapply(temp, function(s) s[2])
temp = strsplit(vars.gnomad$Uploaded_variation, "_")
vars.gnomad$variant = sapply(temp, function(s) s[1])
vars.gnomad$gnomad_freqs = sapply(temp, function(s) s[2])

vars.mgrb = vars.mgrb[,colnames(vars.mgrb) != "Uploaded_variation"]
vars.gnomad = vars.gnomad[,colnames(vars.gnomad) != "Uploaded_variation"]
vars = merge(vars.mgrb, vars.gnomad, all = TRUE)
rm(vars.mgrb, vars.gnomad, temp, infile.mgrb, infile.gnomad)

colnames(vars) = tolower(colnames(vars))

temp = strsplit(vars$variant, ":")
vars$chrom = sapply(temp, function(s) s[1])
vars$pos = as.integer(sapply(temp, function(s) s[2]))
vars$ref = sapply(temp, function(s) s[3])
vars$alt = sapply(temp, function(s) s[4])
temp = strsplit(vars$mgrb_freqs, ":")
vars$AC.MGRB = as.integer(sapply(temp, function(s) s[2]))
vars$AN.MGRB = as.integer(sapply(temp, function(s) s[3]))
vars$NAA.MGRB = as.integer(sapply(temp, function(s) s[4]))
temp = strsplit(vars$gnomad_freqs, ":")
vars$AC.GnomAD = as.integer(sapply(temp, function(s) s[2]))
vars$AN.GnomAD = as.integer(sapply(temp, function(s) s[3]))
vars$NAA.GnomAD = as.integer(sapply(temp, function(s) s[4]))
rm(temp)

# Keep only autosomal variants
vars = vars[grepl("^[0-9]+$", vars$chrom),]

# Impute biallelic genotype frequencies
vars$NS.MGRB = vars$AN.MGRB/2 
vars$NS.GnomAD = vars$AN.GnomAD/2
vars$NRA.MGRB = vars$AC.MGRB - 2*vars$NAA.MGRB
vars$NRA.GnomAD = vars$AC.GnomAD - 2*vars$NAA.GnomAD
vars$NRR.MGRB = vars$NS.MGRB - vars$NRA.MGRB - vars$NAA.MGRB
vars$NRR.GnomAD = vars$NS.GnomAD - vars$NRA.GnomAD - vars$NAA.GnomAD

vars = vars[,c("chrom", "pos", "ref", "alt", "symbol", "feature", "ccds", "flags", "consequence", "exon", "hgvsp", 
    "NS.MGRB", "AC.MGRB", "AN.MGRB", "NRR.MGRB", "NRA.MGRB", "NAA.MGRB", 
    "NS.GnomAD", "AC.GnomAD", "AN.GnomAD", "NRR.GnomAD", "NRA.GnomAD", "NAA.GnomAD")]

stopifnot(any(vars$NRR.GnomAD != as.integer(vars$NRR.GnomAD), na.rm = TRUE) == FALSE)
stopifnot(any(vars$NRA.GnomAD != as.integer(vars$NRA.GnomAD), na.rm = TRUE) == FALSE)
stopifnot(any(vars$NAA.GnomAD != as.integer(vars$NAA.GnomAD), na.rm = TRUE) == FALSE)
stopifnot(any(vars$NRR.MGRB != as.integer(vars$NRR.MGRB), na.rm = TRUE) == FALSE)
stopifnot(any(vars$NRA.MGRB != as.integer(vars$NRA.MGRB), na.rm = TRUE) == FALSE)
stopifnot(any(vars$NAA.MGRB != as.integer(vars$NAA.MGRB), na.rm = TRUE) == FALSE)

# Classify variants according to the following table, where the highest hit takes priority:
# Classification   Consequences
# Nonsense         stop_gained
#                  start_lost
#                  splice_donor_variant
#                  splice_acceptor_variant
# Missense         missense_variant
#                  stop_lost
# Synonymous       synonymous_variant
#                  stop_retained_variant
# Unclassified     coding_sequence_variant
#                  splice_region_variant

vars$is_nonsense = grepl("stop_gained|start_lost|splice_donor_variant|splice_acceptor_variant", vars$consequence)
vars$is_missense = grepl("missense_variant|stop_lost", vars$consequence)
vars$is_syn = (grepl("synonymous_variant|stop_retained_variant", vars$consequence) & !(vars$is_nonsense | vars$is_missense))
vars$is_nonsyn = vars$is_nonsense | vars$is_missense
vars$class = NA
vars$class[vars$is_syn] = "Synonymous"
vars$class[vars$is_nonsyn] = "Nonsynonymous"

vars = vars[!is.na(vars$class),]

vars$in_mgrb = !is.na(vars$AC.MGRB) & vars$AC.MGRB > 0
vars$in_gnomad = !is.na(vars$AC.GnomAD) & vars$AC.GnomAD > 0
vars$cohort = "Neither"
vars$cohort[vars$in_mgrb & !vars$in_gnomad] = "MGRB_only"
vars$cohort[!vars$in_mgrb & vars$in_gnomad] = "GnomAD_only"
vars$cohort[vars$in_mgrb & vars$in_gnomad] = "Both"

table(vars$class, vars$cohort)
#                  Both GnomAD_only MGRB_only
#  Nonsynonymous 138802      549540    145095
#  Synonymous     93581      282076     69736
#  N:S            1.483       1.948     2.081

# Not unsurprising that shared variants are more likely to be synonymous than
# unshared -- these will have larger AFs on average.

# The discrepancy in the cohort-private variants is statistically significant:
fisher.test(table(vars$class[vars$cohort != "Both"], vars$cohort[vars$cohort != "Both"]))
# 
#         Fisher's Exact Test for Count Data
# 
# data:
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.9269241 0.9458895
# sample estimates:
# odds ratio
#  0.9363377

# How to design a test for differential dN/dS that takes into account
# variable denominator across loci?  Permutation test?

# Data to illustrate the problem:
# vars[vars$symbol == "KRAS",c("NS.MGRB", "AC.MGRB", "AN.MGRB", "NRR.MGRB", "NRA.MGRB", "NAA.MGRB", "NS.GnomAD", "AC.GnomAD", "AN.GnomAD", "NRR.GnomAD", "NRA.GnomAD", "NAA.GnomAD", "class")]
#        NS.MGRB AC.MGRB AN.MGRB NRR.MGRB NRA.MGRB NAA.MGRB NS.GnomAD AC.GnomAD AN.GnomAD NRR.GnomAD NRA.GnomAD NAA.GnomAD         class
# 123494      NA      NA      NA       NA       NA       NA      7505         1     15010       7504          1          0 Nonsynonymous
# 123495      NA      NA      NA       NA       NA       NA      7482         1     14964       7481          1          0 Nonsynonymous
# 123496      NA      NA      NA       NA       NA       NA      7493         1     14986       7492          1          0 Nonsynonymous
# 123497    2570       2    5140     2568        2        0      7504         4     15008       7500          4          0 Nonsynonymous
# 123498      NA      NA      NA       NA       NA       NA      7500         1     15000       7499          1          0 Nonsynonymous
# 123499    2570       1    5140     2569        1        0      7503         2     15006       7501          2          0 Nonsynonymous
# 123500      NA      NA      NA       NA       NA       NA      7499         1     14998       7498          1          0    Synonymous
# 123501      NA      NA      NA       NA       NA       NA      7504         1     15008       7503          1          0    Synonymous
# 123502    2570       3    5140     2567        3        0      7506         7     15012       7499          7          0    Synonymous
# 123503    2570       1    5140     2569        1        0        NA        NA        NA         NA         NA         NA    Synonymous
# 123504      NA      NA      NA       NA       NA       NA      7503         2     15006       7501          2          0    Synonymous
# 123505      NA      NA      NA       NA       NA       NA      7489         1     14978       7488          1          0    Synonymous
# 123506      NA      NA      NA       NA       NA       NA      7483         1     14966       7482          1          0    Synonymous
# 123507    2570    5140    5140        0        0     2570      7502     15004     15004          0          0       7502    Synonymous
# 123508      NA      NA      NA       NA       NA       NA      7493         1     14986       7492          1          0    Synonymous

# The challenge is that in cases where a variant was not detected in a cohort (NA blocks),
# we don't know the number of samples genotyped in that cohort (ie NS and AN unknown).
# Consequently, the combined cohort NS or AN needed for a permutation test is also not known.

# Can we reasonably assume that AN.MGRB ~ AN.GnomAD, for shared variants?  If this is the
# case, we can impute missing NS or AN using the values from the cohort in which it was
# present.

# The other (more robust) alternative is to focus on shared variants only.  This will be
# particularly stringent, as these shared variants are more likely to be of high AF, and
# therefore low (or even beneficial) effect.  Power may be too low, especially considering 
# the relative scarcity of these variants.

# Visual inspection of AN.MGRB ~ AN.GnomAD indicates little correlation.  For shared variants, 
# AN.MGRB is nearly always > 5000, whereas AN.GnomAD is far more variable.  Imputation would
# be best performed 'within-cohort', by sampling AN marginally from all ANs in that cohort.

# To support this rather strong choice, restrict variants to those with NS / max(NS) > 0.98
# in both cohorts:
temp.sel = (is.na(vars$NS.MGRB) | vars$NS.MGRB / max(vars$NS.MGRB, na.rm = TRUE) > 0.98) &
    (is.na(vars$NS.GnomAD) | vars$NS.GnomAD / max(vars$NS.GnomAD, na.rm = TRUE) > 0.98)
mean(temp.sel)
# [1] 0.9787071
vars = vars[temp.sel,]

table(vars$class, vars$cohort)
#                  Both GnomAD_only MGRB_only
#  Nonsynonymous 138226      532665    145084
#  Synonymous     93193      272700     69732
#  N:S            1.483       1.953     2.081

# As expected given ~ 2% of variants lost, little change.



source("03_variantburden.R")



vb = VariantBurden(vars$symbol, vars$class, vars$AC.MGRB, vars$AN.MGRB, vars$NAA.MGRB, vars$AC.GnomAD, vars$AN.GnomAD, vars$NAA.GnomAD, vars)


# Load variant lists
universe = unique(vb$data$symbol)
# ACMG SF 2.0
acmg = c("APC", "MYH11", "ACTA2", "TMEM43", "DSP", "PKP2", "DSG2", "DSC2", "BRCA1", "BRCA2", "SCN5A", "RYR2", "LMNA", "MYBPC3", "COL3A1", "GLA", "APOB", 
    "LDLR", "MYH7", "TPM1", "MYBPC3", "PRKAG2", "TNNI3", "MYL3", "MYL2", "ACTC1", "RET", "PCSK9", "BMPR1A", "SMAD4", "TNNT2", "TP53", "TGFBR1", "TGFBR2", 
    "TGFBR1", "TGFBR2", "SMAD3", "KCNQ1", "KCNH2", "SCN5A", "MLH1", "MSH2", "MSH6", "PMS2", "RYR1", "CACNA1S", "FBN1", "TGFBR1", "MEN1", "RET", "MUTYH", 
    "NF2", "OTC", "SDHD", "SDHAF2", "SDHC", "SDHB", "STK11", "MUTYH", "PTEN", "RB1", "TSC1", "TSC2", "VHL", "WT1", "ATP7B")
# COSMIC CGC Germline 20171020
cgc = c("ALK","APC","APOBEC3B","AR","ATM","ATR","AXIN2","BAP1","BLM","BMPR1A","BRCA1","BRCA2","BRIP1","BUB1B","CDC73","CDH1","CDK4","CDKN1B","CDKN2A",
    "CHEK2","CXCR4","CYLD","DDB2","DICER1","EGFR","ERBB4","ERCC2","ERCC3","ERCC4","ERCC5","EXT1","EXT2","FANCA","FANCC","FANCD2","FANCE","FANCF","FANCG",
    "FAT1","FH","FLCN","GPC3","HNF1A","HRAS","KDR","KIT","LMO1","LZTR1","MAX","MEN1","MLH1","MPL","MSH2","MSH6","MUTYH","NBN","NF1","NF2","PALB2","PDGFRA",
    "PHOX2B","PMS2","POLE","PRF1","PRKAR1A","PTCH1","PTEN","PTPN13","RB1","RECQL4","RET","SBDS","SDHA","SDHAF2","SDHB","SDHC","SDHD","SETBP1","SMAD4",
    "SMARCB1","SMARCE1","SPOP","STAT3","STK11","SUFU","TERT","TGFBR2","TMEM127","TP53","TP63","TSC1","TSC2","TSHR","VHL","WAS","WRN","WT1","XPA","XPC")
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4767849/ table S1
cardiac = c("ABCC9", "ABCG8", "ACTA1", "ACTA2", "ACTC1", "ACTN2", "AKAP9", "ANK2", "ANKRD1", "APOB", "APOE", "BAG3", "BRAF", "CACNA1C", "CACNA2D1", 
    "CACNB2", "CALM1", "CALR3", "CASQ2", "CAV3", "CBL", "CETP", "COL3A1", "COX15", "CRYAB", "CSRP3", "CTF1", "DES", "DMD", "DNAJC19", "DOLK", "DSC2", 
    "DSG2", "DSP", "DTNA", "EFEMP2", "ELN", "EMD", "FBN1", "FHL1", "FHL2", "FKTN", "FXN", "GAA", "GATAD1", "GLA", "GPD1L", "HCN4", "ILK", "JPH2", "JUP", 
    "KCND3", "KCNE1", "KCNE2", "KCNE3", "KCNH2", "KCNJ2", "KCNJ5", "KCNJ8", "KCNQ1", "KLF10", "KRAS", "LAMA2", "LAMA4", "LAMP2", "LDB3", "LDLR", "LDLRAP1", 
    "LMNA", "LTBP2", "MAP2K1", "MAP2K2", "MIB1", "CAVIN4", "MYBPC3", "MYH11", "MYH6", "MYH7", "MYL2", "MYL3", "MYLK", "MYLK2", "MYO6", "MYOZ2", "MYPN", 
    "NEXN", "NOTCH1", "NPPA", "NRAS", "PCSK9", "PDLIM3", "PKP2", "PLN", "PRDM16", "PRKAG2", "PTPN11", "RAF1", "RANGRF", "RBM20", "RYR2", "SCN2B", "SCN3B", 
    "SCN4B", "SCN5A", "SCO2", "SGCB", "SGCD", "SHOC2", "SLC25A4", "SLC2A10", "SMAD3", "SNTA1", "SOS1", "SREBF2", "TAZ", "TCAP", "TGFB2", "TGFB3", "TGFBR1", 
    "TGFBR2", "TMEM43", "TMPO", "TNNC1", "TNNI3", "TNNT2", "TPM1", "TRDN", "TRIM63", "TRPM4", "TTN", "TXNRD2", "VCL")
alldisease = unique(c(acmg, cgc, cardiac))
controls = setdiff(universe, alldisease)
set.seed(314159)
randoms = replicate(1000, as.character(sample(universe, length(alldisease))), simplify = FALSE)


# Prepare for multicore computation
set.seed(314159, kind = "L'Ecuyer-CMRG")
library(plyr)
library(doParallel)
registerDoParallel(28)


# Perform bootstrapping of statistics
boot.acmg = boot.VariantBurden(vb, test_symbols = acmg, control_symbols = controls, R = 200, seed = 314159, ncpus = 28, parallel = "multicore")
boot.cgc = boot.VariantBurden(vb, test_symbols = cgc, control_symbols = controls, R = 200, seed = 314159, ncpus = 28, parallel = "multicore")
boot.cardiac = boot.VariantBurden(vb, test_symbols = cardiac, control_symbols = controls, R = 200, seed = 314159, ncpus = 28, parallel = "multicore")
boot.all = boot.VariantBurden(vb, test_symbols = alldisease, control_symbols = controls, R = 200, seed = 314159, ncpus = 28, parallel = "multicore")
boot.randoms.t0 = laply(randoms, function(random) boot.VariantBurden(vb, test_symbols = random, control_symbols = setdiff(universe, random), R = 0, seed = 314159)$t0, .parallel = TRUE)


# Perform permutation tests
test.acmg = test_ns_tc.VariantBurden(vb, test_symbols = acmg, control_symbols = controls, B = 1000, seed = 314159, .parallel = TRUE)
test.cgc = test_ns_tc.VariantBurden(vb, test_symbols = cgc, control_symbols = controls, B = 1000, seed = 314159, .parallel = TRUE)
test.cardiac = test_ns_tc.VariantBurden(vb, test_symbols = cardiac, control_symbols = controls, B = 1000, seed = 314159, .parallel = TRUE)
test.all = test_ns_tc.VariantBurden(vb, test_symbols = alldisease, control_symbols = controls, B = 1000, seed = 314159, .parallel = TRUE)
test.randoms.pval = laply(randoms[1:100], function(random) { pval = test_ns_tc.VariantBurden(vb, test_symbols = random, control_symbols = setdiff(universe, random), B = 1000, seed = 314159)$p; cat("."); flush.console(); pval }, .parallel = TRUE)


# Save results
saveRDS(list(
    vb = vb,
    genelists = list(acmg = acmg, cgc = cgc, cardiac = cardiac, all = alldisease, randoms = randoms),
    boots = list(acmg = boot.acmg, cgc = boot.cgc, cardiac = boot.cardiac, all = boot.all, randoms.t0 = boot.randoms.t0),
    tests = list(acmg = test.acmg, cgc = test.cgc, cardiac = test.cardiac, all = test.all, randoms.pval = test.randoms.pval)
), file = "03_summary.rds")


# Plots
svg("03_summarise_%02d.svg", height = 8, width = 8)
plotboot.VariantBurden(boot.acmg, c("C1TN/C1TS", "C1CN/C1CS", "C2TN/C2TS", "C2CN/C2CS"))
plotboot.VariantBurden(boot.cgc, c("C1TN/C1TS", "C1CN/C1CS", "C2TN/C2TS", "C2CN/C2CS"))
plotboot.VariantBurden(boot.cardiac, c("C1TN/C1TS", "C1CN/C1CS", "C2TN/C2TS", "C2CN/C2CS"))
plotboot.VariantBurden(boot.all, c("C1TN/C1TS", "C1CN/C1CS", "C2TN/C2TS", "C2CN/C2CS"))
hist(log(boot.randoms.t0[,"((C1TN/C1TS)/(C1CN/C1CS)) / ((C2TN/C2TS)/(C2CN/C2CS))"]), breaks = 20)
hist(test.randoms.pval)
plot(ecdf(test.randoms.pval))
dev.off()
