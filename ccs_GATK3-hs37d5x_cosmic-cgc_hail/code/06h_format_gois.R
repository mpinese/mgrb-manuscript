options(java.parameters = "-Xmx4G")

args = commandArgs(trailingOnly = TRUE)
infile = args[[1]]
outfile = args[[2]]

variants = read.table(infile, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gois = read.table("tmp/gois.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

gois$List = gsub(" +", "_", gsub(",", " ", gois$List))

for (class in sort(unique(gois$Class)))
{
    variants = cbind(variants, variants$Symbol %in% gois$Genes[gois$Class == class])
    colnames(variants)[ncol(variants)] = sprintf("Class.%s", class)
}

for (list in sort(unique(gois$List)))
{
    variants = cbind(variants, variants$Symbol %in% gois$Genes[gois$List == list])
    colnames(variants)[ncol(variants)] = sprintf("List.%s", list)
}

library(xlsx)

write.xlsx2(variants, outfile, row.names = FALSE)

