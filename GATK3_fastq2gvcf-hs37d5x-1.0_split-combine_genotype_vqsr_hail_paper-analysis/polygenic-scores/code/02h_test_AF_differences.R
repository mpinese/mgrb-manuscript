# Compare the relative frequency of alleles in MGRB vs GnomAD, split by:
#   * All alleles
#   * All GWAS-significant alleles (EBI GWAS DB)
#   * GWAS-significant alleles not associated with disease
#   * GWAS-significant alleles associated with disease depleted in MGRB

options(echo = TRUE)

pdf("02_AF_differences.pdf")

set.seed(314159)

# Load AFs and GWAS results
afs = readRDS("../MGRB_GnomAD_AFs_outerjoined_common.rds")
ebi_gwas = readRDS("../gwas_catalog_v1.0.1-associations_e91_r2018-03-13.cleaned.rds")
manual_gwas = read.table("../manual_polygenic_scores.GnomAD_NFE_AFs.models", stringsAsFactors = FALSE, header = TRUE, sep = "\t")


afs$isWellGenotyped = afs$rate_MGRB >= 0.98 & afs$rate_gnomad >= 0.98
afs = afs[afs$isWellGenotyped,]


afs$vid = sprintf("%d:%d:%s:%s", afs$chrom, afs$pos, afs$ref, afs$alt)


# Process the EBI GWAS results to extract:
#   * All alleles
#   * All 'negative control' (not associated with disease) alleles

traits.negcon.morphology = c(
    "body height",
    "eye morphology measurement",
    "facial morphology",
    "facial morphology measurement",
    "hair morphology",
    "lip morphology measurement",
    "lip morphology measurement, chin morphology measurement",
    "mouth morphology measurement",
    "nose morphology measurement",
    "nose morphology measurement, mouth morphology measurement")

traits.negcon.psych = c(
    "intelligence, self reported educational attainment",
    "self reported educational attainment",
    "economic and social preference",
    "neuroticism measurement",
    "intelligence, specific language impairment, dyslexia",
    "non-word reading",
    "word reading",
    "episodic memory",
    "word list delayed recall measurement, memory performance")

ebi_gwas.all = ebi_gwas
ebi_gwas.negcon.morphology = ebi_gwas[tolower(ebi_gwas$trait) %in% traits.negcon.morphology,]
ebi_gwas.negcon.psych = ebi_gwas[tolower(ebi_gwas$trait) %in% traits.negcon.psych,]


# Process the manual GWAS entries to extract GWAS-significant alleles linked
# to conditions we expect to be depleted in MGRB.

# Drop unhelpful scores.  These are Wood (expect neutral), and Deelen
# (all coefficients unidirectional, so doesn't fit in test framework).
manual_gwas.MGRB_depleted = manual_gwas[!(manual_gwas$id %in% c("Height:Wood:10.1038/ng.3097", "Longevity:Deelen:10.1093/hmg/ddu139")),]




# Define a test function to compare the allele rates.
depletion_test = function(vids, coefs, afs, seed = 314159, method = c("exact", "approx", "fisher"))
{
    # Null hypothesis: Whether or not an allele is enriched in MGRB vs GnomAD
    # is independent of coef.
    method = match.arg(method)
    if (method == "exact")
        stopifnot(require(Exact))

    saved_seed = .Random.seed
    set.seed(seed)

    # 1. In the case of duplicate vids, randomly select one.
    perm = sample.int(length(vids))
    vids = vids[perm]
    coefs = coefs[perm]
    drop = duplicated(vids)
    vids = vids[!drop]
    coefs = coefs[!drop]

    # 2. Extract AFs for the vids
    message("Extracting AFs...")
    afs_mgrb = afs$AF_MGRB[match(vids, afs$vid)]
    afs_gnomad = afs$AF_gnomad[match(vids, afs$vid)]

    # 3. Drop vids with invalid values
    drop = coefs == 0 | is.na(coefs) | is.na(afs_mgrb) | is.na(afs_gnomad) | afs_mgrb == afs_gnomad
    vids = vids[!drop]
    coefs = coefs[!drop]
    afs_mgrb = afs_mgrb[!drop]
    afs_gnomad = afs_gnomad[!drop]

    # 4. Define variants as enriched in MGRB vs not, and positive coef vs not.
    coef.pos = coefs > 0
    afs.mgrb_enrich = afs_mgrb > afs_gnomad

    # 5. Construct a contingency table 
    #                        MGRB AF
    #                  Depleted  Enriched
    # Coef  Negative
    #       Positive
    tbl = matrix(
        c(sum(!coef.pos & !afs.mgrb_enrich), sum(coef.pos & !afs.mgrb_enrich), sum(!coef.pos & afs.mgrb_enrich), sum(coef.pos & afs.mgrb_enrich)), nrow = 2,
        dimnames = list("Coefficient" = c("Negative", "Positive"), "MGRB" = c("Depleted", "Enriched")))

    message("Testing...")
    if (method == "exact")
        test = exact.test(tbl, alternative = "two.sided", method = "boschloo", model = "Binomial", cond.row = TRUE)
    else if (method == "prop")
        test = prop.test(tbl[,2], tbl[,1] + tbl[,2])
    else
        test = fisher.test(tbl)

    assign(".Random.seed", saved_seed, envir = .GlobalEnv)

    list(table = tbl, test = test)
}


depletion_test(ebi_gwas.negcon.morphology$vid, ebi_gwas.negcon.morphology$logbeta, afs, method = "fisher")
depletion_test(ebi_gwas.negcon.psych$vid, ebi_gwas.negcon.psych$logbeta, afs, method = "fisher")
#depletion_test(ebi_gwas.all$vid, ebi_gwas.all$logbeta, afs, method = "approx")
depletion_test(manual_gwas.MGRB_depleted$vid, manual_gwas.MGRB_depleted$coef, afs, method = "fisher")

# The below is not useful -- almost all tests are underpowered
# library(plyr)
# dlply(manual_gwas.MGRB_depleted, .(id), function(d) depletion_test(d$vid, d$coef, afs), .progress = "text")




enrichment_plot = function(afs, vids = NULL)
{
    library(plyr)
    library(ggplot2)
    if (!is.null(vids))
        afs = afs[afs$vid %in% vids,]
    AF_combined_bin = cut((afs$AF_gnomad + afs$AF_MGRB)/2, breaks = seq(0, 1, 0.05))
    temp.af_enrich = ddply(afs, .(AF_combined_bin), function(d) c(n_enriched_MGRB = sum(d$enriched_MGRB), n = nrow(d)))
    ggplot(temp.af_enrich, aes(x = AF_combined_bin, y = n_enriched_MGRB/n)) + geom_point() + ylim(0, 1) + theme_bw() + ylab("Rate enriched in MGRB") + xlab("Allele frequency") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
}


enrichment_plot(afs) + ggtitle("All variants")
#enrichment_plot(afs, ebi_gwas.all$vid) + ggtitle("EBI GWAS variants")
#enrichment_plot(afs, manual_gwas.MGRB_depleted$vid) + ggtitle("Manual GWAS variants")


# Compare MGRB and GnomAD AFs
set.seed(314159)
temp.sample = sample.int(nrow(afs), 1e6)
plot(AF_MGRB ~ AF_gnomad, afs[temp.sample,], pch = ".", col = rgb(0, 0, 0, 0.05), xlab = "GnomAD NFE VAF", ylab = "MGRB VAF")
abline(0, 1, lty = "dotted", col = "yellow")


afs_ebi = afs[afs$vid %in% ebi_gwas$vid,]
wilcox.test(afs_ebi$AF_MGRB - afs_ebi$AF_gnomad)
#binom.test(sum((afs_ebi$AF_MGRB - afs_ebi$AF_gnomad) > 0), nrow(afs_ebi))

dev.off()

