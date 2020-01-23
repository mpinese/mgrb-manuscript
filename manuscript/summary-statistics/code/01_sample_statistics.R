meta = read.csv("C:/Users/markp/Desktop/Projects/MGRB/misc_data/MGRB_Phase2_metadata.csv", stringsAsFactors = FALSE)
meta = meta[meta$Phase2.SampleTier <= 2 & !is.na(meta$Phase2.SampleTier) & meta$Phase2.RelatednessDropped == 0,]
meandepth = read.table("C:/Users/markp/Desktop/Projects/MGRB/MGRB.phase2.meandepth.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE, col.names = c("sampleID", "meandepth"))

meta$cohort = c("A" = "ASPREE", "B" = "45andUp", "Z" = "CRM")[substr(meta$sampleID, 1, 1)]
meta = merge(meta, meandepth)

table(meta$cohort, meta$isFemale)

library(plyr)
ddply(meta, .(cohort), function(d) fivenum(d$YOB))
ddply(meta, .(cohort), function(d) fivenum(d$SBPMean))
ddply(meta, .(cohort), function(d) fivenum(d$HtMtrs))
ddply(meta, .(cohort), function(d) fivenum(d$WtKgs))
ddply(meta, .(cohort), function(d) fivenum(d$AbdoCircCms))
ddply(meta, .(cohort), function(d) fivenum(d$GlcmmolL))
ddply(meta, .(cohort), function(d) mean(d$treatedForHighBP))
ddply(meta, .(cohort), function(d) mean(d$treatedForHighChol))
ddply(meta, .(cohort), function(d) fivenum(d$meandepth))

