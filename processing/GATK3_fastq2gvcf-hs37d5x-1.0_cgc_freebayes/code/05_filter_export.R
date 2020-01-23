library(data.table)

vars = readRDS("../MGRB_phase2.dupmarked.realigned.recalibrated.cgc.freebayes.split.vtnorm.merged.rds")
vep = fread("tmp/vep.out", skip = "#Uploaded_variation")
cosmic = fread("../04_cosmic_83.tsv")

colnames(vep)[1] = "VID"
vars$VID = sprintf("%s:%d:%s:%s", vars$chrom, vars$pos, vars$ref, vars$alt)
cosmic$VID = sprintf("%s:%d:%s:%s", cosmic$chrom, cosmic$pos, cosmic$ref, cosmic$alt)

vars = merge(vars, vep, by = "VID", all.x = TRUE, all.y = FALSE)
vars = merge(vars, cosmic, by = "VID", all.x = TRUE, all.y = FALSE)
rm(cosmic, vep)

vars$csq.syn = grepl("synonymous_variant|stop_retained_variant|start_retained_variant", vars$Consequence)
vars$csq.nonsyn = grepl("missense_variant|splice_acceptor_variant|splice_donor_variant|start_lost|stop_gained|stop_lost|frameshift_variant|inframe_deletion|inframe_insertion", vars$Consequence)
vars$csq.conflicting = vars$csq.syn & vars$csq.nonsyn

mean(is.na(vars$Location))
sum(vars$csq.conflicting)

vars = vars[!is.na(vars$Location) & (vars$csq.syn | vars$csq.nonsyn) & !vars$csq.conflicting,]

temp = head(vars)

vars[vars == "-"] = NA


# eps = 1e-2
# vars$lprRA = pbinom(vars$ad, vars$ad + vars$rd, 0.5, log = TRUE)
# vars$lprAA = pbinom(vars$ad, vars$ad + vars$rd, 1-eps/3, log = TRUE)
# vars$lprAB = pbinom(vars$ad, vars$ad + vars$ado, 0.5, log = TRUE)
# vars$lprRR = pbinom(vars$ad + vars$ado, vars$ad + vars$ado + vars$rd, eps, log = TRUE)


