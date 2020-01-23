setwd("C:/Users/mpinese/Downloads")

library(openxlsx)

dat.depleted = read.xlsx("Supplementary Data 3.xlsx", sheet=4)
dat.anthro = read.xlsx("Supplementary Data 3.xlsx", sheet=6)
dat.behav = read.xlsx("Supplementary Data 3.xlsx", sheet=8)

# Recode MGRB so it's always after the comparator cohorts, for consistency
dat.depleted$VAF.higher.in.MGRB.vs.gnomAD[dat.depleted$VAF.higher.in.MGRB.vs.gnomAD == "MGRB"] = "Z_MGRB"
dat.depleted$VAF.higher.in.MGRB.vs.UKBB[dat.depleted$VAF.higher.in.MGRB.vs.UKBB == "MGRB"] = "Z_MGRB"
dat.anthro$VAF.higher.in.MGRB.vs.gnomAD[dat.anthro$VAF.higher.in.MGRB.vs.gnomAD == "MGRB"] = "Z_MGRB"
dat.anthro$VAF.higher.in.MGRB.vs.UKBB[dat.anthro$VAF.higher.in.MGRB.vs.UKBB == "MGRB"] = "Z_MGRB"
dat.behav$VAF.higher.in.MGRB.vs.gnomAD[dat.behav$VAF.higher.in.MGRB.vs.gnomAD == "MGRB"] = "Z_MGRB"
dat.behav$VAF.higher.in.MGRB.vs.UKBB[dat.behav$VAF.higher.in.MGRB.vs.UKBB == "MGRB"] = "Z_MGRB"

tbl.depleted.gnomad = table(dat.depleted$Association.with.depleted.phenotype, dat.depleted$VAF.higher.in.MGRB.vs.gnomAD)
tbl.depleted.ukbb = table(dat.depleted$Association.with.depleted.phenotype, dat.depleted$VAF.higher.in.MGRB.vs.UKBB)
tbl.anthro.gnomad = table(dat.anthro$Association.with.anthropometric.phenotype, dat.anthro$VAF.higher.in.MGRB.vs.gnomAD)
tbl.anthro.ukbb = table(dat.anthro$Association.with.anthropometric.phenotype, dat.anthro$VAF.higher.in.MGRB.vs.UKBB)
tbl.behav.gnomad = table(dat.behav$Association.with.behavioural.phenotype, dat.behav$VAF.higher.in.MGRB.vs.gnomAD)
tbl.behav.ukbb = table(dat.behav$Association.with.behavioural.phenotype, dat.behav$VAF.higher.in.MGRB.vs.UKBB)

fisher.test(tbl.depleted.gnomad)
fisher.test(tbl.depleted.ukbb)
fisher.test(tbl.anthro.gnomad)
fisher.test(tbl.anthro.ukbb)
fisher.test(tbl.behav.gnomad)
fisher.test(tbl.behav.ukbb)
