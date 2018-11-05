# A large number of SNPs in data/manual_polygenic_scores.GnomAD_NFE_AFs.models are 
# not reliably genotyped (not in the high confidence regions defined by 
# data/MGRBphase2final_dpge15_98pct.GnomAD201_dpge15_98pct.GiaB332HC.bed.xz).  Try
# and rescue these using tag SNPs which are in the HCR.  Use an R^2 threshold of
# 0.9 as assessed in 1000G Phase 3 NFE populations, with the one exception of
# 19:45411941:T:C (APOE rs429358) -- no tag SNP could be identified at this high
# threshold, so 19:45415713:G:A (R^2 = 0.75191) was manually selected instead.

library(GenomicRanges)

# Load models and high confidence genotyping regions

temp.models = read.table("data/manual_polygenic_scores.GnomAD_NFE_AFs.models", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
temp.hcrs = read.table("data/MGRBphase2final_dpge15_98pct.GnomAD201_dpge15_98pct.GiaB332HC.bed.xz", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

offsets = temp.models[temp.models$vid == "OFFSET",]
temp.models = temp.models[temp.models$vid != "OFFSET",]

vid2coord = function(vids)
{
    spl = strsplit(vids, ":")
    list(
        contig = sapply(spl, function(x) x[1]),
        pos = as.integer(sapply(spl, function(x) x[2])),
        ref = sapply(spl, function(x) x[3]),
        alt = sapply(spl, function(x) x[4]))
}

temp.split = vid2coord(temp.models$vid)
models = GRanges(
    seqnames = Rle(temp.split$contig), 
    ranges = IRanges(start = temp.split$pos, end = temp.split$pos + nchar(temp.split$ref) - 1), 
    strand = "*", 
    vid = temp.models$vid, id = temp.models$id, ref = temp.split$ref, alt = temp.split$alt, coef = temp.models$coef, aaf = temp.models$aaf)

hcrs = GRanges(
    seqnames = Rle(temp.hcrs[,1]),
    ranges = IRanges(start = temp.hcrs[,2] + 1, end = temp.hcrs[,3]),
    strand = "*")

rm(temp.models, temp.hcrs)


# Mark model variants as in or out of the HCRs.
models$in_hcr = countOverlaps(models, hcrs) > 0


# Try and rescue variants that are not in the HCRs.  Do this as follows:
# 1. Generate a 1000 genomes NFE haplotype set for imputation:
#    Download 1000G data
#      cd orig
#      parallel wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi ::: $(seq 1 22)
#      parallel wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ::: $(seq 1 22)
#      wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
#      cd ..
#    Subset to biallelic SNPs in NFE populations
#      awk '($2 == "GBR") { print $1 }' < orig/integrated_call_samples_v3.20130502.ALL.panel > gbr_samples.txt
#      awk '($2 == "CEU") { print $1 }' < orig/integrated_call_samples_v3.20130502.ALL.panel > ceu_samples.txt
#      awk '($2 == "IBS") { print $1 }' < orig/integrated_call_samples_v3.20130502.ALL.panel > ibs_samples.txt
#      awk '($2 == "TSI") { print $1 }' < orig/integrated_call_samples_v3.20130502.ALL.panel > tsi_samples.txt
#      cat gbr_samples.txt ceu_samples.txt ibs_samples.txt tsi_samples.txt > nfe_samples.txt
#      parallel bcftools view -a -c 1 -O z -S nfe_samples.txt -m2 -M2 -v snps -o ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.nfe_biallelic_snps.vcf.gz orig/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ::: $(seq 1 22)
#    Relabel with Hail-style vids, and fix duplicated variants (eg 12:8400000:T:G)
#      parallel gzip -dc ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.nfe_biallelic_snps.vcf.gz '|' awk "'"'BEGIN {FS="\t"; OFS="\t"} /^#/; /^[^#]/ { $3 = $1 ":" $2 ":" $4 ":" $5; if (vids[$3] == 0) { vids[$3] = 1; print } }'"'" '>' ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.nfe_biallelic_snps.nodups.vcf ::: $(seq 1 22)
#    Convert to plink format
#      parallel plink --vcf ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.nfe_biallelic_snps.nodups.vcf --a2-allele ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.nfe_biallelic_snps.nodups.vcf 4 3 "'"'#'"'" --make-bed --out ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.nfe_biallelic_snps ::: $(seq 1 22)
#    Merge
#      rm -f merge_files.txt
#      seq 2 22 | while read chr; do echo ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.nfe_biallelic_snps >> merge_files.txt; done
#      plink --bfile ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.nfe_biallelic_snps --merge-list merge_files.txt --make-bed --out ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.nfe_biallelic_snps
#    Clean up
#      rm ALL.chr*.phase3* merge_files.txt
#    GT rate and HWE filtering
#      plink --bfile ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.nfe_biallelic_snps --geno 0.02 --out filter_geno98 --write-snplist
#      plink --bfile ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.nfe_biallelic_snps --keep <(awk '{print $1, $1}' < gbr_samples.txt) --hwe 0.001 --out filter_hwegbr001 --write-snplist
#      plink --bfile ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.nfe_biallelic_snps --keep <(awk '{print $1, $1}' < ceu_samples.txt) --hwe 0.001 --out filter_hweceu001 --write-snplist
#      plink --bfile ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.nfe_biallelic_snps --keep <(awk '{print $1, $1}' < ibs_samples.txt) --hwe 0.001 --out filter_hweibs001 --write-snplist
#      plink --bfile ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.nfe_biallelic_snps --keep <(awk '{print $1, $1}' < tsi_samples.txt) --hwe 0.001 --out filter_hwetsi001 --write-snplist
#      sort filter_geno98.snplist filter_hwegbr001.snplist filter_hweceu001.snplist filter_hweibs001.snplist filter_hwetsi001.snplist | uniq -c | awk '($1 == 5) { print $2 }' > filter_allpass.snplist
#      plink --bfile ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.nfe_biallelic_snps --extract filter_allpass.snplist --make-bed --out data/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.nfe_biallelic_qcpass_snps
#    Clean up
#      rm filter_geno* filter_hwe* ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.nfe_biallelic_snps.* *_samples.txt
# 2. For each model variant not in the HCRs, query for tag snps in 1000G.
#      plink --bfile ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.nfe_biallelic_qcpass_snps --show-tags target_snps.txt --tag-r2 0.9 --tag-kb 250 --list-all --out target_snps
# 3. Import back into R and filter the tag snps themselves for presence in the HCR.
# 4. For any tag snps that remain, choose the closest in physical location to the target snp, and then detect and correct for phase change.

write.table(sort(unique(models$vid[models$in_hcr == FALSE])), "tmp/target_snps.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

system("plink --bfile data/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.nfe_biallelic_qcpass_snps --show-tags tmp/target_snps.txt --tag-r2 0.9 --tag-kb 250 --list-all --out tmp/target_snps")

tags = read.table("tmp/target_snps.tags.list", header = TRUE, stringsAsFactors = FALSE)
tags$TAGS[tags$TAGS %in% c("", "NONE")] = NA

# Find the best tag for each target SNP
find_best_tag = function(tags, target)
{
    if (is.na(tags))
        return(NA)
    coords = vid2coord(strsplit(tags, "\\|")[[1]])
    coords_ranges = GRanges(seqnames = Rle(coords$contig), ranges = IRanges(start = coords$pos, width = nchar(coords$ref)), strand = "*", ref = coords$ref, alt = coords$alt)
    coords_ranges_hcr = coords_ranges[countOverlaps(coords_ranges, hcrs) > 0]
    if (length(coords_ranges_hcr) == 0)
        return(NA)
    target = vid2coord(target)
    target_range = GRanges(seqnames = Rle(target$contig), ranges = IRanges(start = target$pos, width = nchar(target$ref)), strand = "*")
    dists = distance(target_range, coords_ranges_hcr, select = "all")
    best_tag = coords_ranges_hcr[which.min(dists)]
    sprintf("%s:%d:%s:%s", seqnames(best_tag), start(best_tag), best_tag$ref, best_tag$alt)
}
tags$best_tag = mapply(find_best_tag, tags$TAGS, tags$SNP)
tags = tags[,c("SNP", "best_tag")]
colnames(tags)[1] = "target"

# Now determine whether the tag SNP is in phase with the target or not.
# If it's in phase, no changes need to be made to the model.
# If it's out of phase, the model beta must be inverted.
# OFFSETs should be fixed also, but as we're really only interested in
# relative scores leave the OFFSET correction out for now.

query_phase = function(tag, target)
{
    if (is.na(tag))
        return(NA)
    system(sprintf("plink --bfile data/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.nfe_biallelic_qcpass_snps --ld %s %s | grep -Eo 'In phase alleles are [ACGT]{2}/[ACGT]{2}' | sed 's/In phase alleles are //'", tag, target), intern = TRUE)
}

library(plyr)
library(doParallel)
registerDoParallel(28)
tags = ddply(tags, .(target), function(d) {
    message(d$target)
    temp = try(query_phase(d$best_tag, d$target))
    if (class(temp) == "try-error" || length(temp) > 1)     # length > 1 can occur in the case of multiple solns, eg 12:56647911:G:C 12:56647243:C:T ==> GC/CT or GT/CC
        d$phase = NA
    else
        d$phase = temp
    d }, .parallel = TRUE)

# Phase now contains tagA-targetA/tagB-targetB

tags$in_phase = daply(tags, .(target), function(d) {
    if (is.na(d$phase))
        return(NA)
    temp = vid2coord(d$target)
    targ.ref = temp$ref
    targ.alt = temp$alt
    temp = vid2coord(d$best_tag)
    tag.ref = temp$ref
    tag.alt = temp$alt

    inphase1 = sprintf("%s%s/%s%s", tag.ref, targ.ref, tag.alt, targ.alt)
    inphase2 = sprintf("%s%s/%s%s", tag.alt, targ.alt, tag.ref, targ.ref)
    outphase1 = sprintf("%s%s/%s%s", tag.ref, targ.alt, tag.alt, targ.ref)
    outphase2 = sprintf("%s%s/%s%s", tag.alt, targ.ref, tag.ref, targ.alt)

    stopifnot(d$phase %in% c(inphase1, inphase2, outphase1, outphase2))
    d$phase %in% c(inphase1, inphase2)
})


# Map this information back to the models

models$tag_substituted = models$vid %in% tags$target[!is.na(tags$in_phase)]
models$original_vid = NA
models$original_coef = NA
models$original_aaf = NA
models$tag_in_phase = NA
models$original_vid[models$tag_substituted] = models$vid[models$tag_substituted]
models$original_coef[models$tag_substituted] = models$coef[models$tag_substituted]
models$original_aaf[models$tag_substituted] = models$aaf[models$tag_substituted]
models$tag_in_phase[models$tag_substituted] = tags$in_phase[match(models$vid[models$tag_substituted], tags$target)]

models$vid[models$tag_substituted] = tags$best_tag[match(models$vid[models$tag_substituted], tags$target)]
models$coef[models$tag_substituted & models$tag_in_phase == FALSE] = -models$coef[models$tag_substituted & models$tag_in_phase == FALSE]
models$aaf[models$tag_substituted] = NA   # We don't bother extracting the gnomAD AFs here as they won't be used

# Drop variants which have not been rescued by this process
length(models)                                  # 2449
sum(models$in_hcr)                              # 1876
sum(models$in_hcr | models$tag_substituted)     # 2222
models = models[models$in_hcr | models$tag_substituted]

# Save the results
newmodels = data.frame(
    id = models$id,
    vid = models$vid,
    coef = models$coef,
    tag_substituted = models$tag_substituted,
    tag_in_phase = models$tag_in_phase,
    original_vid = models$original_vid,
    original_coef = models$original_coef)

write.csv(newmodels, file = "data/manual_polygenic_scores.hcr_tag_rescued.csv", row.names = FALSE, quote = FALSE)

