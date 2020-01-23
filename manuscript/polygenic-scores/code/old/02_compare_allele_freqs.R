# Get all variants of interest and their betas
models = read.table("../manual_polygenic_scores.GnomAD_NFE_AFs.models", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
models = models[models$vid != "OFFSET",]

# Decompose the variant IDs
temp = strsplit(models$vid, ":")
models$chrom = sapply(temp, function(x) as.integer(x[1]))
models$pos = sapply(temp, function(x) as.integer(x[2]))
models$ref = sapply(temp, function(x) x[3])
models$alt = sapply(temp, function(x) x[4])

# Ensure all variants are autosomal
stopifnot(any(is.na(models$chrom)) == FALSE)
stopifnot(models$chrom >= 1 & models$chrom <= 22)


# Load MGRB genotypes and calculate AFs
mgrb.GT = read.table("../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.popstar.modelonly.dosages.tab", stringsAsFactors = FALSE, header = TRUE, sep = "\t")
mgrb.GT = mgrb.GT[mgrb.GT[,1] %in% models$vid,]
mgrb.vid = mgrb.GT[,1]

mgrb.GT = as.matrix(mgrb.GT[,-1])
mgrb.AC = rowSums(mgrb.GT, na.rm = TRUE)
mgrb.AN = rowSums(!is.na(mgrb.GT))*2      # Safe to assume autosomal due to check above

mgrb.AFs = data.frame(vid = mgrb.vid, mgrb.AC = mgrb.AC, mgrb.AN = mgrb.AN)

rm(mgrb.GT, mgrb.vid, mgrb.AC, mgrb.AN)


# Load GnomAD AFs
library(DBI)
db = dbConnect(RSQLite::SQLite(), "../gnomad.genomes.r2.0.1.sites.combined.split.minrep.NFE.dbSNP.autosomes.dba")

query = dbSendQuery(db, 'SELECT * FROM dbsnp WHERE chrom = :chrom AND pos = :pos AND ref = :ref AND alt = :alt')
dbBind(query, param = list(chrom = as.character(models$chrom), pos = models$pos, ref = models$ref, alt = models$alt))
gnomad.data = dbFetch(query)
dbClearResult(query)
dbDisconnect(db)

gnomad.AFs = data.frame(vid = sprintf("%s:%d:%s:%s", gnomad.data$chrom, gnomad.data$pos, gnomad.data$ref, gnomad.data$alt), gnomad.AC = gnomad.data$AC, gnomad.AN = gnomad.data$AN)

rm(gnomad.data, db, query)


# Merge all data
data = merge(merge(models, mgrb.AFs, by = "vid", all = TRUE), gnomad.AFs, by = "vid", all = TRUE)
rm(models, mgrb.AFs, gnomad.AFs)


# Restrict to variants well genotyped in both cohorts.
# Numeric thresholds on AN correspond to 90% of the cohorts genotyped.
data.complete = data[!is.na(data$mgrb.AN) & !is.na(data$gnomad.AN),]
data.complete = data.complete[data.complete$mgrb.AN >= 4630 & data.complete$gnomad.AN >= 13514,]

nrow(data)
# [1] 1238
nrow(data.complete)
# [1] 1110

temp = data[!(data$vid %in% data.complete$vid),]
temp[order(abs(temp$coef)),]

data.orig = data
data = data.complete
rm(data.complete)



# Test for AF shift in MGRB vs GnomAD
# Classify AF shifts in MGRB - GnomAD as either risk-increasing, risk-decreasing, or neutral.
# Base this on classifying models as either risk-positive, risk-negative, or neutral.
unique(data$id)
#  [1] "Height:Wood:10.1038/ng.3097"                          "BreastCancer:Michailidou:10.1038/nature24284"         "EOCAD:Theriault:10.1161/circgen.117.001849"
#  [4] "ProstateCancer:Hoffmann:10.1158/2159-8290.CD-15-0315" "SystolicBP:Warren:10.1038/ng.3768"                    "DiastolicBP:Warren:10.1038/ng.3768"
#  [7] "BasalCellCarcinoma:Chahal:10.1038/ncomms12510"        "ColorectalCancer:Schumacher:10.1038/ncomms8138"       "AlzheimersDisease:Lambert:10.1038/ng.2802"
# [10] "PulsePressure:Warren:10.1038/ng.3768"                 "Melanoma:Law:10.1038/ng.3373"                         "SquamousCellCarcinoma:Chahal:10.1038/ncomms12048"
# [13] "Longevity:Deelen:10.1093/hmg/ddu139"

# Almost all models are risk-positive: higher score => worse health
data$risk_association = +1
# Height is probably close to neutral
data$risk_association[data$id == "Height:Wood:10.1038/ng.3097"] = 0
# The Deelen Longevity model is risk-negative: higher score => longer life
data$risk_association[data$id == "Longevity:Deelen:10.1093/hmg/ddu139"] = -1

# Classify AF shifts as MGRB higher (positive) or GnomAD higher (negative)
data$af_shift = sign((data$mgrb.AC / data$mgrb.AN) - (data$gnomad.AC / data$gnomad.AN))

# Classify alleles as either risk-associated (positive) or risk-depleting (negative)
data$allele_association = data$risk_association * sign(data$coef)

# The null hypothesis: there should be no af_shift sign bias between the allele_association classes.
# In other words, the probability of an allele being more prevalent in MGRB vs GnomAD should be
# independent of its association with risk.

# Exclude the neutral alleles as they don't contribute much to this comparison
data.non_neutral = data[data$allele_association != 0,]
data.neutral = data[data$allele_association == 0,]
table(data.non_neutral$af_shift, data.non_neutral$allele_association)

# Looks interesting.  The correct test to use here is paired binomial comparison,
# as one margin (allele_association) is fixed, but the other (af_shift) is free.
temp.mgrb_more_common.protective = sum(data.non_neutral$af_shift > 0 & data.non_neutral$allele_association < 0, na.rm = TRUE)
temp.mgrb_more_common.deleterious = sum(data.non_neutral$af_shift > 0 & data.non_neutral$allele_association > 0, na.rm = TRUE)
temp.gnomad_more_common.protective = sum(data.non_neutral$af_shift < 0 & data.non_neutral$allele_association < 0, na.rm = TRUE)
temp.gnomad_more_common.deleterious = sum(data.non_neutral$af_shift < 0 & data.non_neutral$allele_association > 0, na.rm = TRUE)
temp.total.protective = temp.mgrb_more_common.protective + temp.gnomad_more_common.protective
temp.total.deleterious = temp.mgrb_more_common.deleterious + temp.gnomad_more_common.deleterious
prop.test(c(temp.mgrb_more_common.protective, temp.mgrb_more_common.deleterious), c(temp.total.protective, temp.total.deleterious))

# As a control, test the height PRS.  This is large enough to have decent power,
# and we expect no preference for height +ve vs -ve alleles in MGRB vs GnomAD.
temp.mgrb_more_common.tall = sum(data.neutral$af_shift > 0 & data.neutral$coef > 0, na.rm = TRUE)
temp.mgrb_more_common.short = sum(data.neutral$af_shift > 0 & data.neutral$coef < 0, na.rm = TRUE)
temp.gnomad_more_common.tall = sum(data.neutral$af_shift < 0 & data.neutral$coef > 0, na.rm = TRUE)
temp.gnomad_more_common.short = sum(data.neutral$af_shift < 0 & data.neutral$coef < 0, na.rm = TRUE)
temp.total.tall = temp.mgrb_more_common.tall + temp.gnomad_more_common.tall
temp.total.short = temp.mgrb_more_common.short + temp.gnomad_more_common.short
prop.test(c(temp.mgrb_more_common.tall, temp.mgrb_more_common.short), c(temp.total.tall, temp.total.short))

# Hmm... borderline significant for MGRB being depleted for tall variants.
# Although the effect is far more subtle than for the disease phenotypes, could there be a source of bias?

# In all cases it seems that positive effects are associated with lower AF in MGRB.

