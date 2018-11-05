#!/usr/bin/env Rscript
options(echo = TRUE)
args = commandArgs(trailingOnly = TRUE)
message(args[[1]])

# Get the full list of complete pass 2 output files
temp.pipe = pipe(sprintf("find '%s/pass2/shards' -name '*.shard*.done'", args[[1]]))
temp = scan(temp.pipe, what = character(0))
close(temp.pipe)
temp = strsplit(temp, "/")
infiles = data.frame(sample = sapply(temp, function(x) x[2]))
temp = strsplit(sapply(temp, function(x) x[5]), "\\.")
infiles$shard = gsub("^shard", "", sapply(temp, function(x) x[4]))
infiles$chrom = gsub("^chrom", "", sapply(temp, function(x) x[5]))
infiles$start1 = as.integer(gsub("-.*", "", gsub("^pos", "", sapply(temp, function(x) x[6]))))
infiles$end1 = as.integer(gsub(".*-", "", gsub("^pos", "", sapply(temp, function(x) x[6]))))
infiles$variant_file = sprintf("../%s/pass2/shards/%s.soma-snv.pass2.shard%s.chrom%s.pos%d-%d.variants.tsv", infiles$sample, infiles$sample, infiles$shard, infiles$chrom, infiles$start1, infiles$end1)
infiles$background_file = sprintf("../%s/pass2/shards/%s.soma-snv.pass2.shard%s.chrom%s.pos%d-%d.background.tsv", infiles$sample, infiles$sample, infiles$shard, infiles$chrom, infiles$start1, infiles$end1)

sample = infiles$sample[[1]]

suppressPackageStartupMessages(library(plyr))

variants = ddply(infiles, .(shard), function(d)
    read.table(
        file = d$variant_file, 
        colClasses = c(character(0), integer(0), character(0), character(0), integer(0), integer(0), integer(0), integer(0), numeric(0), numeric(0)),
        col.names = c("chrom", "pos", "context", "alt", "strand", "dp", "rd", "ad", "sens", "snr"), 
        stringsAsFactors = FALSE
    )
)
variants$alt = as.character(variants$alt)
variants$alt[variants$alt == "TRUE"] = "T"

if (nrow(variants) > 0) {
    variants$sample = sample
} else {
    variants$sample = character(0)
}

backgrounds = ddply(infiles, .(shard), function(d)
    read.table(
        file = d$background_file, 
        colClasses = c(character(0), numeric(0)),
        col.names = c("context", "background"), 
        stringsAsFactors = FALSE
    )
)
backgrounds$sample = sample
backgrounds = backgrounds[,c("sample", "context", "background")]

saveRDS(variants, sprintf("%s/pass2/%s.combined.variants.rds", args[[1]], sample))
saveRDS(backgrounds, sprintf("%s/pass2/%s.combined.backgrounds.rds", args[[1]], sample))

# Compute normalised signature burden
signatures = expand.grid(ref_upstream = c("A", "C", "G", "T"), ref = c("C", "T"), ref_downstream = c("A", "C", "G", "T"), alt = c("A", "C", "G", "T"))
signatures = signatures[as.character(signatures$ref) != signatures$alt,]
signatures$signature = sprintf("%s%s %s.%s", signatures$ref, signatures$alt, signatures$ref_upstream, signatures$ref_downstream)
signatures$context = paste(signatures$ref_upstream, signatures$ref, signatures$ref_downstream, sep = "")
signatures$count = unlist(alply(signatures, 1, function(d) sum(variants$context == d$context & variants$alt == d$alt)))
signatures$background = unlist(alply(signatures, 1, function(d) sum(backgrounds$background[backgrounds$context == d$context])))
signatures$burden = signatures$count / signatures$background
signatures$sample = sample
signatures = signatures[,c("sample", "context", "alt", "signature", "count", "background", "burden")]

saveRDS(signatures, sprintf("%s/pass2/%s.combined.burden.rds", args[[1]], sample))
