#!/usr/bin/env Rscript

AUTHOR = "Mark Pinese <m.pinese@garvan.org.au>"
VERSION = "0.0.3 (26 Oct 2017)"

suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(RSQLite))
suppressPackageStartupMessages(library(DBI))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(dbplyr))


read_population_associations = function(assoc_file, ancestry_file, pop_regex = "european")
{
    ancestry = read.table(ancestry_file, sep = "\t", header = TRUE, comment.char = "", quote = "", row.names = NULL, stringsAsFactors = FALSE)
    colnames(ancestry) = c(colnames(ancestry)[-1], "dummy")
    ancestry = ancestry[,-ncol(ancestry)]

    count.popn.initial = daply(ancestry, .(STUDY.ACCCESSION), function(a) sum(a$NUMBER.OF.INDIVDUALS[a$STAGE == "initial" & grepl(pop_regex, tolower(a$BROAD.ANCESTRAL.CATEGORY))]))
    count.popn.replication = daply(ancestry, .(STUDY.ACCCESSION), function(a) sum(a$NUMBER.OF.INDIVDUALS[a$STAGE == "replication" & grepl(pop_regex, tolower(a$BROAD.ANCESTRAL.CATEGORY))]))

    assoc = read.table(assoc_file, sep = "\t", header = TRUE, comment.char = "", quote = "", stringsAsFactors = FALSE)
    assoc = assoc[,c("PUBMEDID", "STRONGEST.SNP.RISK.ALLELE", "PVALUE_MLOG", "P.VALUE..TEXT.", "OR.or.BETA", "MAPPED_TRAIT", "MAPPED_TRAIT_URI", "STUDY.ACCESSION", "RISK.ALLELE.FREQUENCY")]
    assoc$pop = pop_regex
    assoc$pop.n.initial = count.popn.initial[assoc$STUDY.ACCESSION]
    assoc$pop.n.replication = count.popn.replication[assoc$STUDY.ACCESSION]

    # Calculate ln(Beta), discarding invalid values
    assoc$log_beta = log(assoc$OR.or.BETA)
    assoc = assoc[!(is.na(assoc$log_beta)) & is.finite(assoc$log_beta),,drop=FALSE]

    assoc$RISK.ALLELE.FREQUENCY = as.numeric(assoc$RISK.ALLELE.FREQUENCY)

    # Discard non-SNP associations, or rows with incomplete data (missing at
    # least one of rsid, effect allele, or beta)
    temp = strsplit(assoc$STRONGEST.SNP.RISK.ALLELE, "-")
    assoc$rsid = sapply(temp, function(l) l[1])
    assoc$effect_allele = toupper(sapply(temp, function(l) l[2]))
    assoc = assoc[!is.na(assoc$effect_allele) & !is.na(assoc$OR.or.BETA) & nchar(assoc$rsid) > 2 & substr(assoc$rsid, 1, 2) == "rs" & grepl("^[ACGT]$", assoc$effect_allele),]
    assoc = assoc[,colnames(assoc) != "STRONGEST.SNP.RISK.ALLELE"]

    assoc$id = paste(assoc$STUDY.ACCESSION, gsub(":", ";", assoc$MAPPED_TRAIT), gsub("^ *\\(", "", gsub(" *\\)$", "", assoc$P.VALUE..TEXT)), sep = ":")

    assoc = assoc[,c("id", "STUDY.ACCESSION", "PUBMEDID", "MAPPED_TRAIT", "MAPPED_TRAIT_URI", "pop", "pop.n.initial", "pop.n.replication", "rsid", "effect_allele", "RISK.ALLELE.FREQUENCY", "log_beta", "PVALUE_MLOG")]
    colnames(assoc) = c(    "id", "accession",       "pubmed",   "trait",        "trait_uri",        "pop", "n.initial",     "n.replication",     "rsid", "effect_allele", "effect_allele_freq", "log_beta", "log10_p")

    assoc[!duplicated(assoc),]
}


get_snps_from_db = function(rsids, dbsnp_tbl)
{
    # Fetch dbSNP data for the rsids in assoc from a dbsnp-AF database

    snps.matching_rsids = as.data.frame(select(filter(rsids, name %in% rsids), 
        chrom, chromEnd, name, refNCBI, strand, observed, alleles))
    colnames(snps.matching_rsids) = c("chrom", "pos", "ref", "alt", "rsid", "AC", "AN")
    message(sprintf("%d rsIDs found in database", length(unique(snps.matching_rsids$rsid))))

    # Filter down to biallelic unambiguous strand SNPs.
    snps.matching_rsids = snps.matching_rsids[snps.matching_rsids$observed %in% c("A/C", "A/G", "C/T", "G/T"),,drop=FALSE]
    message(sprintf("%d strand-specific SNPs", length(unique(snps.matching_rsids$rsid))))

    # Filter down to canonical chromosomes
    snps.matching_rsids$chrom = sub("^chr", "", snps.matching_rsids$chrom)
    snps.matching_rsids = snps.matching_rsids[grepl("^([0-9]+|X|Y)$", snps.matching_rsids$chrom),,drop=FALSE]
    message(sprintf("%d rsIDs in canonical chromosomes", length(unique(snps.matching_rsids$rsid))))

    # Extract alleles
    temp = t(simplify2array(strsplit(snps.matching_rsids$alleles, ",")))
    colnames(temp) = c("allele1", "allele2")
    snps.matching_rsids = cbind(snps.matching_rsids, temp)
    snps.matching_rsids$allele1 = as.character(snps.matching_rsids$allele1)
    snps.matching_rsids$allele2 = as.character(snps.matching_rsids$allele2)

    snps.matching_rsids
}


is_valid_strand_specific_snp = function(allele1, allele2)
{
    (allele1 != allele2) && 
    (allele1 %in% c("A", "C", "G", "T")) &&
    (allele2 %in% c("A", "C", "G", "T")) &&
    (
        (allele1 == "A" && allele2 != "T") ||
        (allele1 == "T" && allele2 != "A") ||
        (allele1 == "C" && allele2 != "G") ||
        (allele1 == "G" && allele2 != "C")
    )
}


check_snp_against_db = function(target_rsid, effect_allele, dbsnp_tbl)
{
    dbsnp_matches = dbsnp_tbl %>% filter(rsid == target_rsid)

    # Run a gauntlet of consistency checks.

    # Fail on missing SNPs or multi-allelics
    if (as.data.frame(count(dbsnp_matches))$n != 1)
        return(list(passes = FALSE, chrom = NA, pos = NA, ref = NA, alt = NA, aaf = NA))

    dbsnp_matches = as.data.frame(dbsnp_matches)

    # Fail on ambiguous or invalid SNPs.
    if (!is_valid_strand_specific_snp(dbsnp_matches$ref, dbsnp_matches$alt))
        return(list(passes = FALSE, chrom = NA, pos = NA, ref = NA, alt = NA, aaf = NA))

    # Fail on SNPs which don't have the effect allele
    if (effect_allele != dbsnp_matches$ref && effect_allele != dbsnp_matches$alt)
        return(list(passes = FALSE, chrom = NA, pos = NA, ref = NA, alt = NA, aaf = NA))

    return(list(
        passes = TRUE,
        chrom = dbsnp_matches$chrom, pos = dbsnp_matches$pos, 
        ref = dbsnp_matches$ref, alt = dbsnp_matches$alt, 
        aaf = dbsnp_matches$AC / dbsnp_matches$AN))
}


rebase_association = function(target_rsid, effect_allele, effect_allele_freq, log_beta, dbsnp_tbl)
{
    # Run the SNP through a series of checks, and retrieve information on the SNP.
    dbsnp_check_result = check_snp_against_db(target_rsid, effect_allele, dbsnp_tbl)

    if (dbsnp_check_result$passes == TRUE)
    {
        # SNP passes checks.  Rebase to the reference if needed
        if (effect_allele == dbsnp_check_result$ref)
        {
            # effect allele is the ref: rebasing required.
            # The idea behind rebasing is that we have a predictor of the following form:
            #    y = a1*x + a0
            # and wish to express it in terms of x' = (2-x).  This necessitates a change
            # to a1 and a0:
            #    a1*x + a0 = a1'*(2-x) + a0'  for all x
            # => a1*x + a0 - 2*a1' + x*a1' - a0' = 0
            # => a1 = -a1'
            # => a1*x + a0 + 2*a1 - a1*x - a0' = 0
            # => a0 - a0' + 2*a1 = 0
            # => a0' = 2*a1 + a0
            # (Note the above is only valid for autosomes, due to the 2-x).
            # (Further note that in this case a0 = 0, as the EBI GWAS database does not
            # store offset terms.)
            rebased.log_beta = -log_beta
            rebased.offset = 2*log_beta
        }
        else
        {
            # effect allele is the alt: no rebasing required
            rebased.log_beta = log_beta
            rebased.offset = 0
        }
    }
    else
    {
        # Failure.  We add a term to the offset to marginalise the missing SNP (this
        # is a rough imputation).  Using the notation above:
        #   y = a1*x + a0
        # But we do not have x.  We substitute xhat = 2*EAF, where EAF is the effect
        # allele frequency in the derivation cohort.  Again this is autosome-specific.
        # Then, y = a1*xhat + a0 = a1*2*EAF + a0.
        # We include the a1*2*EAF into a modified a0:
        #   a0' = a0 + a1*2*EAF
        # so that the missing variant is now rolled into the offset term.  Again a0
        # is initially zero.
        rebased.log_beta = 0
        rebased.offset = log_beta*2*effect_allele_freq
    }

    # Update the result structure with the rebased values and return
    dbsnp_check_result$log_beta = rebased.log_beta
    dbsnp_check_result$offset = rebased.offset
    return(dbsnp_check_result)
}


rebase_associations = function(assoc, dbsnp_tbl)
{
    message(sprintf("%d rsIDs in input", length(unique(assoc$rsid))))

    rebased_assoc = ddply(assoc, colnames(assoc), function(d) 
        as.data.frame(rebase_association(d$rsid, d$effect_allele, d$effect_allele_freq, d$log_beta, dbsnp_tbl)),
        .progress = "text")

    offsets = dlply(rebased_assoc, .(id), function(d) {
        offset_sum = sum(d$offset, na.rm = TRUE)
        if (is.na(offset_sum))
            offset_sum = 0
        offset_sum
    })
    rebased_assoc = rebased_assoc[rebased_assoc$passes,]

    list(rebased = rebased_assoc[,c("id", "rsid", "chrom", "pos", "ref", "alt", "aaf", "log_beta")], offsets = offsets)
}


drop_small_studies = function(assoc, min.n = 5000)
{
    selected_ids = unique(assoc$id[assoc$n.initial >= min.n & assoc$n.replication >= min.n])

    assoc[assoc$id %in% selected_ids,,drop=FALSE]
}


drop_small_models = function(assoc, min.loci = 5)
{
    loci_per_id = daply(assoc, .(id), nrow)
    selected_ids = names(loci_per_id)[loci_per_id >= min.loci]

    assoc[assoc$id %in% selected_ids,,drop=FALSE]
}


export_rebased_scores = function(assoc.rebased, output_path)
{
    offsets = assoc.rebased$offsets
    coefficients = assoc.rebased$rebased
    coefficients$vid = paste(coefficients$chrom, coefficients$pos, coefficients$ref, coefficients$alt, sep = ":")

    offset_df = data.frame(id = names(offsets), rsid = NA, chrom = NA, pos = NA, ref = NA, alt = NA, aaf = NA, log_beta = unlist(offsets), vid = "OFFSET")
    offset_df = offset_df[offset_df$id %in% unique(coefficients$id),]

    augmented_coefficients = rbind(coefficients, offset_df)
    augmented_coefficients = augmented_coefficients[order(augmented_coefficients$id, augmented_coefficients$chrom, augmented_coefficients$pos),]
    rownames(augmented_coefficients) = NULL

    write.table(augmented_coefficients[,c("id", "vid", "rsid", "aaf", "log_beta")], output_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}



library(docopt)

sprintf('Convert EBI GWAS summary tables to POPSTAR PTS format.

Usage:
  ebi2pts.R [options] <associations> <ancestry> <database>
  ebi2pts.R -h

Options: 
  -n <MIN_N>      Minimum GWAS sample and replication cohort size [default: 5000]
  -m <MIN_M>      Minimum number of loci in each PTS [default: 10]
  -b <POPN>       Target population regular expression [default: european]
  -o <OUT>        Path to output PTS file [default: /dev/stdout]
  -h              Print this help

Arguments: 
  associations  Associations TSV from EBI GWAS Catalog (gwas_catalog_v1.0.1-associations_*)
  ancestry      Ancestry TSV from EBI GWAS Catalog (gwas_catalog-ancestry_*)
  database      dbSNP and allele frequency SQLite database (gnomad.genomes.r2.0.1.sites.combined.split.minrep.NFE.dbSNP.autosomes.dba)

%s
%s
', AUTHOR, VERSION) -> doc

opts = docopt(doc)

if (opts$h)
{
    print(doc)
    stop(0)
}

dbsnp_db_conn = DBI::dbConnect(RSQLite::SQLite(), dbname = opts$database, flags = SQLITE_RO)
dbsnp_tbl = tbl(dbsnp_db_conn, "dbsnp")

associations.raw = read_population_associations(opts$associations, opts$ancestry, opts$b)
message(sprintf("Loaded %d polygenic scores, %d loci total", length(unique(associations.raw$id)), length(unique(associations.raw$rsid))))
if (as.integer(opts$n) > 0)
{
    message("Filtering by cohort size... ", appendLF = FALSE)
    associations.raw = drop_small_studies(associations.raw, min.n = as.integer(opts$n))
    message(sprintf("%d polygenic scores, %d loci total remain", length(unique(associations.raw$id)), length(unique(associations.raw$rsid))))
}
if (as.integer(opts$m) > 1)
{
    message("Filtering by score size... ", appendLF = FALSE)
    associations.raw = drop_small_models(associations.raw, min.loci = as.integer(opts$m))
    message(sprintf("%d polygenic scores, %d loci total remain", length(unique(associations.raw$id)), length(unique(associations.raw$rsid))))
}

message("Rebasing associations... ")
associations.rebased = rebase_associations(associations.raw, dbsnp_tbl)
if (as.integer(opts$m) > 1)
    associations.rebased$rebased = drop_small_models(associations.rebased$rebased, min.loci = as.integer(opts$m))
message(sprintf("Done.  Final count: %d polygenic scores, %d loci total", length(unique(associations.rebased$rebased$id)), length(unique(associations.rebased$rebased$rsid))))
export_rebased_scores(associations.rebased, opts$o)
