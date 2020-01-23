#!/usr/bin/env Rscript
options(echo = TRUE)

# Define the variant signatures
signatures = expand.grid(ref_upstream = c("A", "C", "G", "T"), ref = c("C", "T"), ref_downstream = c("A", "C", "G", "T"), alt = c("A", "C", "G", "T"))
signatures = signatures[as.character(signatures$ref) != signatures$alt,]
signatures = signatures[order(signatures$ref, signatures$alt, signatures$ref_upstream, signatures$ref_downstream),]
signatures$signature = sprintf("%s%s %s.%s", signatures$ref, signatures$alt, signatures$ref_upstream, signatures$ref_downstream)
stopifnot(nrow(signatures) == 96)

# Get the full list of complete pass 2 output files
infiles = sort(scan(pipe("find .. -wholename '*/pass2/*.combined.burden.rds'"), what = character(0)))

# Create signature x sample count and background matrices to populate
samples = gsub("\\..*", "", gsub(".*/", "", infiles))
count = matrix(NA, nrow = nrow(signatures), ncol = length(infiles), dimnames = list(signature = signatures$signature, sample = samples))
background = matrix(NA, nrow = nrow(signatures), ncol = length(infiles), dimnames = list(signature = signatures$signature, sample = samples))

# Populate matrices
for (i in 1:length(infiles))
{
    message(samples[[i]])
    i.data = readRDS(infiles[[i]])
    stopifnot(i.data$sample == samples[[i]])
    stopifnot(i.data$sample == colnames(count)[i])
    stopifnot(i.data$sample == colnames(background)[i])
    count[,i] = i.data$count[match(rownames(count), i.data$signature)]
    background[,i] = i.data$background[match(rownames(background), i.data$signature)]
}

saveRDS(list(count = count, background = background), "../03_mgrb_soma-snv_pass2.combined.burden.rds")
