# Load the EBI GWAS catalog
gwas = read.table("../gwas_catalog_v1.0.1-associations_e91_r2018-03-13.tsv", sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE, comment.char = "")

# Process to strand-specific SNPs with minimal information
gwas = gwas[grepl("rs[0-9]+-[ACGT]", gwas$STRONGEST.SNP.RISK.ALLELE),]
gwas = gwas[grepl("^[0-9]+$", gwas$CHR_ID),]
gwas = gwas[grepl("^[0-9]+$", gwas$CHR_POS),]
temp = strsplit(gwas$STRONGEST.SNP.RISK.ALLELE, "-")
gwas$SNP = sapply(temp, function(x) x[1])
gwas$RISK.ALLELE = sapply(temp, function(x) x[2])
gwas = gwas[gwas$CNV == "N",]
gwas = gwas[,c("STUDY.ACCESSION", "MAPPED_TRAIT", "MAPPED_TRAIT_URI", "CHR_ID", "SNP", "RISK.ALLELE", "PVALUE_MLOG", "OR.or.BETA")]
gwas$MAPPED_TRAIT_URI = gsub(".*/", "", gwas$MAPPED_TRAIT_URI)
gwas$CHR_ID = as.integer(gwas$CHR_ID)
gwas = gwas[!is.na(gwas$OR.or.BETA) & gwas$OR.or.BETA > 0 & gwas$OR.or.BETA != 1,]
gwas = gwas[gwas$MAPPED_TRAIT != "",]
gwas = gwas[grepl("^[ACGT]$", gwas$RISK.ALLELE),]


# Add SNP allele information for biallelic strand-specific SNPs only.
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
gwas$DBSNP.CHROM = NA
gwas$DBSNP.POS = NA
gwas$DBSNP.ALLELE1 = NA
gwas$DBSNP.ALLELE2 = NA
for (temp.chrom in unique(gwas$CHR_ID))
{
    message(temp.chrom)
    temp.snps = as.data.frame(snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, as.character(temp.chrom)))
    temp.snps = temp.snps[temp.snps$alleles_as_ambig %in% c("R", "Y", "K", "M"),]
    temp.sel = gwas$CHR_ID == temp.chrom
    temp.match = match(gwas$SNP[temp.sel], temp.snps$RefSNP_id)
    gwas$DBSNP.CHROM[temp.sel] = as.integer(temp.snps$seqnames[temp.match])
    gwas$DBSNP.POS[temp.sel] = as.integer(temp.snps$pos[temp.match])
    temp.alleles = temp.snps$alleles_as_ambig[temp.match]
    gwas$DBSNP.ALLELE1[temp.sel] = c("R" = "A", "Y" = "C", "K" = "G", "M" = "A")[temp.alleles]
    gwas$DBSNP.ALLELE2[temp.sel] = c("R" = "G", "Y" = "T", "K" = "T", "M" = "C")[temp.alleles]
}

# Remove any SNPs without dbSNP data, or chromosome disagreement between the EBI database and dbSNP.
gwas = gwas[!is.na(gwas$DBSNP.CHROM) & gwas$DBSNP.CHROM == gwas$CHR_ID,]


# Add reference sequence information
library(BSgenome.Hsapiens.UCSC.hg19)
gwas = gwas[order(gwas$DBSNP.CHROM, gwas$DBSNP.POS),]
gwas$REF.ALLELE = NA
for (temp.chrom in unique(gwas$CHR_ID))
{
    message(temp.chrom)
    temp.seq = BSgenome.Hsapiens.UCSC.hg19[[sprintf("chr%d", temp.chrom)]]
    temp.sel = gwas$DBSNP.CHROM == temp.chrom
    gwas$REF.ALLELE[temp.sel] = sapply(gwas$DBSNP.POS[temp.sel], function(pos) { result = try(as.character(temp.seq[pos])); if (class(result) == "try-error") { result = NA }; result })
}

# Remove any loci with bad reference data
gwas = gwas[gwas$REF.ALLELE != "N" & !is.na(gwas$REF.ALLELE),]

# Keep only entries with the reference allele in {allele1, allele2}
gwas = gwas[gwas$REF.ALLELE == gwas$DBSNP.ALLELE1 | gwas$REF.ALLELE == gwas$DBSNP.ALLELE2,]

# Keep only entries with the risk allele in {allele1, allele2}
gwas = gwas[gwas$RISK.ALLELE == gwas$DBSNP.ALLELE1 | gwas$RISK.ALLELE == gwas$DBSNP.ALLELE2,]

# Define the alt allele based on the DBSNP alleles and the ref allele
gwas$ALT.ALLELE = NA
gwas$ALT.ALLELE[gwas$REF.ALLELE == gwas$DBSNP.ALLELE1] = gwas$DBSNP.ALLELE2[gwas$REF.ALLELE == gwas$DBSNP.ALLELE1]
gwas$ALT.ALLELE[gwas$REF.ALLELE == gwas$DBSNP.ALLELE2] = gwas$DBSNP.ALLELE1[gwas$REF.ALLELE == gwas$DBSNP.ALLELE2]
gwas = gwas[!is.na(gwas$ALT.ALLELE),]

# We now have enough data to transform the GWAS data into chrom:pos:ref:alt beta form.
gwas$vid = sprintf("%d:%d:%s:%s", gwas$DBSNP.CHROM, gwas$DBSNP.POS, gwas$REF.ALLELE, gwas$ALT.ALLELE)

# Transform beta to vid orientation
gwas$logbeta = NA
gwas$logbeta[gwas$ALT.ALLELE == gwas$RISK.ALLELE] = log(gwas$OR.or.BETA[gwas$ALT.ALLELE == gwas$RISK.ALLELE])
gwas$logbeta[gwas$REF.ALLELE == gwas$RISK.ALLELE] = -log(gwas$OR.or.BETA[gwas$REF.ALLELE == gwas$RISK.ALLELE])

# Remove unnecessary fields
gwas = gwas[,c("STUDY.ACCESSION", "MAPPED_TRAIT", "MAPPED_TRAIT_URI", "vid", "logbeta")]
colnames(gwas) = c("sid", "trait", "tid", "vid", "logbeta")
gwas = gwas[order(gwas$sid),]


# Save
saveRDS(gwas, "../gwas_catalog_v1.0.1-associations_e91_r2018-03-13.cleaned.rds")
