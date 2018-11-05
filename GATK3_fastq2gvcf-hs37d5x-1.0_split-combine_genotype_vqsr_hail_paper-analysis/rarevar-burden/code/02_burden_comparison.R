library(Biostrings)

# Load and merge variants and annotation
variants = read.csv("../MGRB_GnomAD_SweGen_AFs_forrvas_analysis.csv", 
    col.names = c(
        "contig", "pos", "ref", "alt",
        "nRR_gnomad", "nRA_gnomad", "nAA_gnomad", "nmissing_gnomad",
        "nRR_mgrb_orig", "nRA_mgrb_orig", "nAA_mgrb_orig", "nmissing_mgrb_orig",
        "nRR_mgrb_somfilt", "nRA_mgrb_somfilt", "nAA_mgrb_somfilt", "nmissing_mgrb_somfilt",
        "nRR_swegen", "nRA_swegen", "nAA_swegen", "nmissing_swegen"), 
    stringsAsFactors = FALSE)
annotation = read.table("../MGRB_GnomAD_SweGen_AFs_forrvas_analysis.vep.filtered.tsv", sep = "\t", quote = "", header = TRUE, stringsAsFactors = FALSE)
variants = merge(variants, annotation, by = c("contig", "pos", "ref", "alt"), all = FALSE)
rm(annotation)

# Integrate protein length information to detect far C-term truncations.
temp.ensembl_protein_seqs = readAAStringSet("Homo_sapiens.GRCh38.pep.all.fa.gz")
temp.ensembl_protein_size = data.frame(id = gsub(" .*", "", names(temp.ensembl_protein_seqs)), size = width(temp.ensembl_protein_seqs))
rm(temp.ensembl_protein_seqs)
variants$protein_id = gsub(":.*", "", variants$HGVSp)
variants$protein_length = temp.ensembl_protein_size$size[match(variants$protein_id, temp.ensembl_protein_size$id)]
rm(temp.ensembl_protein_size)

# Filter variants to remove low-confidence cases
variants = variants[!grepl("cds_end_NF|cds_start_NF", variants$FLAGS),]

# Rearrange fields
variants$variant_exon = as.integer(gsub("/.*", "", variants$EXON))
variants$total_exons = as.integer(gsub("[0-9]+/", "", variants$EXON))
variants = variants[,c(
    "contig", "pos", "ref", "alt", "SYMBOL", "Consequence", "Protein_position", "protein_length", "variant_exon", "total_exons", "Feature", "CCDS", "HGVSp", 
    "nRR_gnomad", "nRA_gnomad", "nAA_gnomad", "nmissing_gnomad", 
    "nRR_mgrb_orig", "nRA_mgrb_orig", "nAA_mgrb_orig", "nmissing_mgrb_orig",
    "nRR_mgrb_somfilt", "nRA_mgrb_somfilt", "nAA_mgrb_somfilt", "nmissing_mgrb_somfilt",
    "nRR_swegen", "nRA_swegen", "nAA_swegen", "nmissing_swegen")]
colnames(variants) = c(
    "chrom", "pos", "ref", "alt", "symbol", "consequence", "variant_residue", "total_residues", "variant_exon", "total_exons", "vep_feature", "CCDS", "HGVSp", 
    "nRR_gnomad", "nRA_gnomad", "nAA_gnomad", "nmissing_gnomad", 
    "nRR_mgrb_orig", "nRA_mgrb_orig", "nAA_mgrb_orig", "nmissing_mgrb_orig",
    "nRR_mgrb_somfilt", "nRA_mgrb_somfilt", "nAA_mgrb_somfilt", "nmissing_mgrb_somfilt",
    "nRR_swegen", "nRA_swegen", "nAA_swegen", "nmissing_swegen")
variants$variant_residue = as.integer(gsub("-.*", "", gsub("^\\?-", "", variants$variant_residue)))

# Add a variant class field
variants$variant_class = NA
variants$variant_class[nchar(variants$ref) == 1 & nchar(variants$alt) == 1] = "SNV"
variants$variant_class[nchar(variants$ref) == 1 & nchar(variants$alt) > 1 & nchar(variants$alt) < 20] = "small_insertion"
variants$variant_class[nchar(variants$ref) == 1 & nchar(variants$alt) > 1 & nchar(variants$alt) >= 20] = "large_insertion"
variants$variant_class[nchar(variants$ref) > 1 & nchar(variants$ref) < 20 & nchar(variants$alt) == 1] = "small_deletion"
variants$variant_class[nchar(variants$ref) > 1 & nchar(variants$ref) >= 20 & nchar(variants$alt) == 1] = "large_deletion"
variants$variant_class[nchar(variants$ref) > 1 & nchar(variants$alt) > 1] = "insdel"
variants$variant_class = factor(variants$variant_class, levels = c("SNV", "small_insertion", "large_insertion", "small_deletion", "large_deletion", "insdel"))

# Add a summary effect field
variants$effect = NA
variants$effect[is.na(variants$effect) & grepl("frameshift_variant|stop_gained", variants$consequence) & variants$variant_residue / variants$total_residues < 0.95] = "nonsense"
variants$effect[is.na(variants$effect) & grepl("frameshift_variant|stop_gained", variants$consequence) & variants$variant_residue / variants$total_residues >= 0.95] = "nonsense_terminal"
variants$effect[is.na(variants$effect) & grepl("splice_acceptor_variant|splice_donor_variant", variants$consequence) & variants$variant_residue / variants$total_residues < 0.95] = "splice"
variants$effect[is.na(variants$effect) & grepl("splice_acceptor_variant|splice_donor_variant", variants$consequence) & variants$variant_residue / variants$total_residues >= 0.95] = "splice_terminal"
variants$effect[is.na(variants$effect) & grepl("inframe_deletion|inframe_insertion|protein_altering_variant", variants$consequence)] = "missense_large"
variants$effect[is.na(variants$effect) & grepl("missense_variant", variants$consequence)] = "missense"
variants$effect[is.na(variants$effect) & grepl("synonymous_variant", variants$consequence)] = "synonymous"
variants$effect = factor(variants$effect, levels = c("nonsense", "nonsense_terminal", "splice", "splice_terminal", "missense_large", "missense", "synonymous"))

variants$effect_nosplicing = NA
variants$effect_nosplicing[is.na(variants$effect_nosplicing) & grepl("frameshift_variant|stop_gained", variants$consequence) & variants$variant_residue / variants$total_residues < 0.95] = "nonsense"
variants$effect_nosplicing[is.na(variants$effect_nosplicing) & grepl("frameshift_variant|stop_gained", variants$consequence) & variants$variant_residue / variants$total_residues >= 0.95] = "nonsense_terminal"
variants$effect_nosplicing[is.na(variants$effect_nosplicing) & grepl("inframe_deletion|inframe_insertion|protein_altering_variant", variants$consequence)] = "missense_large"
variants$effect_nosplicing[is.na(variants$effect_nosplicing) & grepl("missense_variant", variants$consequence)] = "missense"
variants$effect_nosplicing[is.na(variants$effect_nosplicing) & grepl("synonymous_variant", variants$consequence)] = "synonymous"
variants$effect_nosplicing = factor(variants$effect_nosplicing, levels = c("nonsense", "nonsense_terminal", "missense_large", "missense", "synonymous"))

# Consequences not covered in the above schema:
# start_lost
# stop_lost


# Collapse down to two syn vs nonsyn calls:
variants$trunc_vs_syn = NA
variants$trunc_vs_syn[variants$effect %in% c("nonsense")] = "truncating"
variants$trunc_vs_syn[variants$effect %in% c("synonymous")] = "synonymous"
variants$nonsyn_vs_syn = NA
variants$nonsyn_vs_syn[variants$effect %in% c("nonsense", "nonsense_terminal", "missense_large", "missense")] = "nonsynonymous"
variants$nonsyn_vs_syn[variants$effect %in% c("synonymous")] = "synonymous"



# Import gene lists.
variants$genelist_ACMG_ADSD = variants$symbol %in% c(
    "ACTA2", "ACTC1", "APC", "APOB", "BMPR1A", "BRCA1", "BRCA2", "CACNA1S", "COL3A1", "DSC2", "DSG2", "DSP", "FBN1", "KCNH2", 
    "KCNQ1", "LDLR", "LMNA", "MEN1", "MLH1", "MSH2", "MSH6", "MYBPC3", "MYH11", "MYH7", "MYL2", "MYL3", "NF2", "OTC", "PCSK9", 
    "PKP2", "PMS2", "PRKAG2", "PTEN", "RB1", "RET", "RYR1", "RYR2", "SCN5A", "SDHAF2", "SDHB", "SDHC", "SDHD", "SMAD3", "SMAD4", 
    "STK11", "TGFBR1", "TGFBR2", "TMEM43", "TNNI3", "TNNT2", "TP53", "TPM1", "TSC1", "TSC2", "VHL", "WT1")
variants$genelist_ACMG_AR = variants$symbol %in% c("ATP7B", "MUTYH")
variants$genelist_ACMG_XL = variants$symbol %in% c("GLA", "OTC")
variants$genelist_CGC_somatic = variants$symbol %in% c(
    "A1CF", "ABI1", "ABL1", "ABL2", "ACKR3", "ACSL3", "ACSL6", "ACVR1", "ACVR2A", "AFF1", "AFF3", "AFF4", "AKAP9", "AKT1", 
    "AKT2", "AKT3", "ALDH2", "ALK", "AMER1", "ANK1", "APC", "AR", "ARAF", "ARHGAP26", "ARHGAP5", "ARHGEF10", "ARHGEF10L", 
    "ARHGEF12", "ARID1A", "ARID1B", "ARID2", "ARNT", "ASPSCR1", "ASXL1", "ASXL2", "ATF1", "ATIC", "ATM", "ATP1A1", "ATP2B3", 
    "ATR", "ATRX", "AXIN1", "AXIN2", "B2M", "BAP1", "BARD1", "BAX", "BAZ1A", "BCL10", "BCL11A", "BCL11B", "BCL2", "BCL2L12",
    "BCL3", "BCL6", "BCL7A", "BCL9", "BCL9L", "BCLAF1", "BCOR", "BCORL1", "BCR", "BIRC3", "BIRC6", "BMP5", "BRAF", "BRCA1", 
    "BRCA2", "BRD3", "BRD4", "BTG1", "BTK", "C15orf65", "C2orf44", "CACNA1D", "CALR", "CAMTA1", "CANT1", "CARD11", "CARS", 
    "CASC5", "CASP3", "CASP8", "CASP9", "CBFA2T3", "CBFB", "CBL", "CBLB", "CBLC", "CCDC6", "CCNB1IP1", "CCNC", "CCND1", 
    "CCND2", "CCND3", "CCNE1", "CCR4", "CCR7", "CD209", "CD274", "CD28", "CD74", "CD79A", "CD79B", "CDC73", "CDH1", "CDH10", 
    "CDH11", "CDH17", "CDK12", "CDK6", "CDKN1A", "CDKN1B", "CDKN2A", "CDKN2C", "CDX2", "CEBPA", "CEP89", "CHCHD7", "CHD2", 
    "CHD4", "CHIC2", "CHST11", "CIC", "CIITA", "CLIP1", "CLP1", "CLTC", "CLTCL1", "CNBD1", "CNBP", "CNOT3", "CNTNAP2", 
    "CNTRL", "COL1A1", "COL2A1", "COL3A1", "COX6C", "CPEB3", "CREB1", "CREB3L1", "CREB3L2", "CREBBP", "CRLF2", "CRNKL1", 
    "CRTC1", "CRTC3", "CSF1R", "CSF3R", "CSMD3", "CTCF", "CTNNA2", "CTNNB1", "CTNND1", "CTNND2", "CUL3", "CUX1", "CXCR4", 
    "CYLD", "CYP2C8", "CYSLTR2", "DAXX", "DCAF12L2", "DCC", "DCTN1", "DDIT3", "DDR2", "DDX10", "DDX3X", "DDX5", "DDX6", 
    "DEK", "DGCR8", "DICER1", "DNAJB1", "DNM2", "DNMT3A", "DROSHA", "DUX4L1", "EBF1", "ECT2L", "EED", "EGFR", "EIF1AX", 
    "EIF3E", "EIF4A2", "ELF3", "ELF4", "ELK4", "ELL", "ELN", "EML4", "EP300", "EPAS1", "EPHA3", "EPHA7", "EPS15", "ERBB2", 
    "ERBB3", "ERBB4", "ERC1", "ERG", "ESR1", "ETNK1", "ETV1", "ETV4", "ETV5", "ETV6", "EWSR1", "EZH2", "EZR", "FAM131B", 
    "FAM135B", "FAM46C", "FAM47C", "FAS", "FAT1", "FAT3", "FAT4", "FBLN2", "FBXO11", "FBXW7", "FCGR2B", "FCRL4", "FES", 
    "FEV", "FGFR1", "FGFR1OP", "FGFR2", "FGFR3", "FGFR4", "FHIT", "FIP1L1", "FKBP9", "FLI1", "FLNA", "FLT3", "FLT4", 
    "FNBP1", "FOXA1", "FOXL2", "FOXO1", "FOXO3", "FOXO4", "FOXP1", "FOXR1", "FSTL3", "FUBP1", "FUS", "GAS7", "GATA1", 
    "GATA2", "GATA3", "GLI1", "GMPS", "GNA11", "GNAQ", "GNAS", "GOLGA5", "GOPC", "GPC5", "GPHN", "GRIN2A", "GRM3", 
    "H3F3A", "H3F3B", "HERPUD1", "HEY1", "HIF1A", "HIP1", "HIST1H3B", "HIST1H4I", "HLA-A", "HLF", "HMGA1", "HMGA2", 
    "HMGN2P46", "HNF1A", "HNRNPA2B1", "HOOK3", "HOXA11", "HOXA13", "HOXA9", "HOXC11", "HOXC13", "HOXD11", "HOXD13", "HRAS", 
    "HSP90AA1", "HSP90AB1", "ID3", "IDH1", "IDH2", "IGF2BP2", "IGH", "IGK", "IGL", "IKBKB", "IKZF1", "IL2", "IL21R", 
    "IL6ST", "IL7R", "IRF4", "IRS4", "ISX", "ITGAV", "ITK", "JAK1", "JAK2", "JAK3", "JAZF1", "JUN", "KAT6A", "KAT6B", 
    "KAT7", "KCNJ5", "KDM5A", "KDM5C", "KDM6A", "KDR", "KDSR", "KEAP1", "KIAA1549", "KIAA1598", "KIF5B", "KIT", "KLF4", 
    "KLF6", "KLK2", "KMT2A", "KMT2C", "KMT2D", "KNSTRN", "KRAS", "KTN1", "LARP4B", "LASP1", "LCK", "LCP1", "LEF1", 
    "LEPROTL1", "LHFP", "LIFR", "LMNA", "LMO1", "LMO2", "LPP", "LRIG3", "LRP1B", "LSM14A", "LYL1", "LZTR1", "MAF", "MAFB", 
    "MALAT1", "MALT1", "MAML2", "MAP2K1", "MAP2K2", "MAP2K4", "MAP3K1", "MAP3K13", "MAPK1", "MAX", "MB21D2", "MDM2", "MDM4", 
    "MDS2", "MECOM", "MED12", "MEN1", "MET", "MGMT", "MITF", "MKL1", "MLF1", "MLH1", "MLLT1", "MLLT10", "MLLT11", "MLLT3", 
    "MLLT4", "MLLT6", "MN1", "MNX1", "MPL", "MSH2", "MSH6", "MSI2", "MSN", "MTCP1", "MTOR", "MUC1", "MUC16", "MUC4", "MYB", 
    "MYC", "MYCL", "MYCN", "MYD88", "MYH11", "MYH9", "MYO5A", "MYOD1", "N4BP2", "NAB2", "NACA", "NBEA", "NCKIPSD", "NCOA1", 
    "NCOA2", "NCOA4", "NCOR1", "NCOR2", "NDRG1", "NF1", "NF2", "NFATC2", "NFE2L2", "NFIB", "NFKB2", "NFKBIE", "NIN", 
    "NKX2-1", "NONO", "NOTCH1", "NOTCH2", "NPM1", "NR4A3", "NRAS", "NRG1", "NSD1", "NT5C2", "NTRK1", "NTRK3", "NUMA1", 
    "NUP214", "NUP98", "NUTM1", "NUTM2A", "NUTM2B", "OLIG2", "OMD", "P2RY8", "PABPC1", "PAFAH1B2", "PAX3", "PAX5", "PAX7", 
    "PAX8", "PBRM1", "PBX1", "PCBP1", "PCM1", "PDCD1LG2", "PDE4DIP", "PDGFB", "PDGFRA", "PDGFRB", "PER1", "PHF6", "PHOX2B", 
    "PICALM", "PIK3CA", "PIK3CB", "PIK3R1", "PIM1", "PLAG1", "PLCG1", "PML", "POLD1", "POLE", "POLG", "POLQ", "POT1", 
    "POU2AF1", "POU5F1", "PPARG", "PPFIBP1", "PPM1D", "PPP2R1A", "PPP6C", "PRCC", "PRDM1", "PRDM16", "PRDM2", "PREX2", 
    "PRKACA", "PRKAR1A", "PRKCB", "PRPF40B", "PRRX1", "PSIP1", "PTCH1", "PTEN", "PTK6", "PTPN11", "PTPN13", "PTPN6", 
    "PTPRB", "PTPRC", "PTPRD", "PTPRK", "PTPRT", "PWWP2A", "QKI", "RABEP1", "RAC1", "RAD17", "RAD21", "RAD51B", "RAF1", 
    "RALGDS", "RANBP2", "RAP1GDS1", "RARA", "RB1", "RBM10", "RBM15", "REL", "RET", "RFWD3", "RGPD3", "RGS7", "RHOA", 
    "RHOH", "RMI2", "RNF213", "RNF43", "ROBO2", "ROS1", "RPL10", "RPL22", "RPL5", "RPN1", "RSPO2", "RSPO3", "RUNDC2A", 
    "RUNX1", "RUNX1T1", "S100A7", "SALL4", "SDC4", "SDHA", "Sep-05", "Sep-06", "Sep-09", "SET", "SETBP1", "SETD1B", 
    "SETD2", "SF3B1", "SFPQ", "SFRP4", "SGK1", "SH2B3", "SH3GL1", "SIRPA", "SIX1", "SIX2", "SKI", "SLC34A2", "SLC45A3", 
    "SMAD2", "SMAD3", "SMAD4", "SMARCA4", "SMARCB1", "SMARCD1", "SMC1A", "SMO", "SND1", "SOCS1", "SOX2", "SOX21", "SPECC1", 
    "SPEN", "SPOP", "SRC", "SRGAP3", "SRSF2", "SRSF3", "SS18", "SS18L1", "SSX1", "SSX2", "SSX4", "STAG1", "STAG2", "STAT3", 
    "STAT5B", "STAT6", "STIL", "STK11", "STRN", "SUFU", "SUZ12", "SYK", "TAF15", "TAL1", "TAL2", "TBL1XR1", "TBX3", "TCEA1", 
    "TCF12", "TCF3", "TCF7L2", "TCL1A", "TEC", "TERT", "TET1", "TET2", "TFE3", "TFEB", "TFG", "TFPT", "TFRC", "TGFBR2", 
    "THRAP3", "TLX1", "TLX3", "TMPRSS2", "TNC", "TNFAIP3", "TNFRSF14", "TNFRSF17", "TOP1", "TP53", "TP63", "TPM3", "TPM4", 
    "TPR", "TRA", "TRAF7", "TRB", "TRD", "TRIM24", "TRIM27", "TRIM33", "TRIP11", "TRRAP", "TSC1", "TSC2", "TSHR", "U2AF1", 
    "UBR5", "USP44", "USP6", "USP8", "VAV1", "VHL", "VTI1A", "WHSC1", "WHSC1L1", "WIF1", "WNK2", "WT1", "WWTR1", "XPO1", 
    "YWHAE", "ZBTB16", "ZCCHC8", "ZEB1", "ZFHX3", "ZMYM3", "ZNF198", "ZNF278", "ZNF331", "ZNF384", "ZNF429", "ZNF479", 
    "ZNF521", "ZNRF3", "ZRSR2")
variants$genelist_ACMG_ADSD_not_somatic = variants$genelist_ACMG_ADSD & !variants$genelist_CGC_somatic
variants$genelist_ACMG_AR_not_somatic = variants$genelist_ACMG_AR & !variants$genelist_CGC_somatic


# Define variants as "rare" if they are present at MAF < 1% in ALL of GnomAD, MGRB, and SweGen
variants$maf_gnomad = (variants$nRA_gnomad + 2*variants$nAA_gnomad) / (2*(variants$nRR_gnomad + variants$nRA_gnomad + variants$nAA_gnomad))
variants$maf_mgrb_orig = (variants$nRA_mgrb_orig + 2*variants$nAA_mgrb_orig) / (2*(variants$nRR_mgrb_orig + variants$nRA_mgrb_orig + variants$nAA_mgrb_orig))
variants$maf_mgrb_somfilt = (variants$nRA_mgrb_somfilt + 2*variants$nAA_mgrb_somfilt) / (2*(variants$nRR_mgrb_somfilt + variants$nRA_mgrb_somfilt + variants$nAA_mgrb_somfilt))
variants$maf_swegen = (variants$nRA_swegen + 2*variants$nAA_swegen) / (2*(variants$nRR_swegen + variants$nRA_swegen + variants$nAA_swegen))
variants$rare_gnomad = (variants$maf_gnomad < 0.01) | is.na(variants$maf_gnomad)
variants$rare_mgrb_orig = (variants$maf_mgrb_orig < 0.01) | is.na(variants$maf_mgrb_orig)
variants$rare_mgrb_somfilt = (variants$maf_mgrb_somfilt < 0.01) | is.na(variants$maf_mgrb_somfilt)
variants$rare_swegen = (variants$maf_swegen < 0.01) | is.na(variants$maf_swegen)
variants$rare_all = variants$rare_gnomad & variants$rare_mgrb_orig & variants$rare_swegen


# Now test the hypotheses.
# We expect that there will be more individuals with a "pathogenic genotype"
# in GnomAD than in MGRB.  Here "pathogenic genotype" is a truncation in
# an ACMG gene, at het or hom alt for AD/SD genes, and hom alt only for AR.
# Do not consider XL.

# Attempt to correct for inevitable platform differences by a dN/dS-like
# normalisation.  For the AD/SD case, normalise by the total of het and
# hom alt synonymous variants; for the AR case, normalise by the count
# of hom alt syn variants only.

# Define our statistics as:
#   S_ADSD = (sum((AA_MGRB_i+RA_MGRB_i)*D_i*ACMG_ADSD_i)/sum((RR_MGRB_i+RA_MGRB_i)*S_i*ACMG_ADSD_i)) / (sum((RR_GnomAD_i+RA_GnomAD_i)*D_i*ACMG_ADSD_i)/sum((RR_GnomAD_i+RA_GnomAD_i)*S_i*ACMG_ADSD_i))
#   S_AR   = (sum(AA_MGRB_i*D_i*ACMG_AR_i)/sum(AA_MGRB_i*S_i*ACMG_AR_i)) / (sum(AA_GnomAD_i*D_i*ACMG_AR_i)/sum(AA_GnomAD_i*S_i*ACMG_AR_i))
# where:
#   AA_MGRB_i   is the frequency of hom alt variants in MGRB at locus i, likewise
#               for RA, GnomAD.
#   D_i         is 1 if variant i is deleterious or pathogenic, else 0.
#   S_i         is 1 if variant i is silent, else 0.
#   ACMG_ADSD_i is 1 if the gene of variant i is an AD or SD ACMG gene, likewise
#               for AR.

# Get approximate bounds on these by bootstrapping.  Bootstrap at both
# the sample level and variant level -- sample for obvious reasons, variant
# as we're only interrogating a subset of loci in the genes of interest.

# As sample-level data are not available, some assumptions need to be made
# to enable sample-level bootstrapping.  Specifically, we assume that all 
# variants are independent, and therefore that resampling variants 
# independently is valid.  As all variants are rare, this is reasonable.


sum(variants$genelist_ACMG_ADSD & variants$trunc_vs_syn == "truncating" & variants$rare_all, na.rm = TRUE)
sum(variants$genelist_ACMG_ADSD & variants$trunc_vs_syn == "synonymous" & variants$rare_all, na.rm = TRUE)
sum(variants$genelist_ACMG_AR & variants$trunc_vs_syn == "truncating" & variants$rare_all, na.rm = TRUE)
sum(variants$genelist_ACMG_AR & variants$trunc_vs_syn == "synonymous" & variants$rare_all, na.rm = TRUE)
sum(variants$genelist_ACMG_ADSD_not_somatic & variants$trunc_vs_syn == "truncating" & variants$rare_all, na.rm = TRUE)
sum(variants$genelist_ACMG_ADSD_not_somatic & variants$trunc_vs_syn == "synonymous" & variants$rare_all, na.rm = TRUE)
sum(variants$genelist_ACMG_AR_not_somatic & variants$trunc_vs_syn == "truncating" & variants$rare_all, na.rm = TRUE)
sum(variants$genelist_ACMG_AR_not_somatic & variants$trunc_vs_syn == "synonymous" & variants$rare_all, na.rm = TRUE)


multi4_sample = function(n1, n2, n3, n4)
{
    # Multinomial resampling by hierarchical binomial draws
    n = n1 + n2 + n3 + n4
    m = length(n)

    f1_1234 = n1/n
    f2_234 = n2/(n-n1)
    f3_34 = n3/(n3+n4)

    f1_1234[n == 0] = 0
    f2_234[n-n1 == 0] = 0
    f3_34[n3+n4 == 0] = 0
    
    r1 = rbinom(m, n, f1_1234)
    r2 = rbinom(m, n - r1, f2_234)
    r3 = rbinom(m, n - r1 - r2, f3_34)
    r4 = as.integer(n - r1 - r2 - r3)

    stopifnot(r4 >= 0)

    list(r1 = r1, r2 = r2, r3 = r3, r4 = r4)
}


redraw_samples = function(nRR, nRA, nAA, nmiss)
{
    na = is.na(nRR) | is.na(nRA) | is.na(nAA) | is.na(nmiss)
    m = length(na)
    
    nRR_resamp = rep(NA, m)
    nRA_resamp = rep(NA, m)
    nAA_resamp = rep(NA, m)
    nmiss_resamp = rep(NA, m)

    resampled = multi4_sample(nRR[!na], nRA[!na], nAA[!na], nmiss[!na])

    nRR_resamp[!na] = resampled$r1
    nRA_resamp[!na] = resampled$r2
    nAA_resamp[!na] = resampled$r3
    nmiss_resamp[!na] = resampled$r4

    list(RR = nRR_resamp, RA = nRA_resamp, AA = nAA_resamp, miss = nmiss_resamp)
}





enrichment_statistic = function(variants, loci, deleterious, test_cohort = "mgrb_orig", ref_cohort = "gnomad", consider_AA = TRUE, consider_RA = TRUE, boot_variants = FALSE, boot_samples = FALSE)
{
    # variants:    data frame as per top level namespace
    # loci:        boolean vector of length equal to rows of variants; TRUE if locus should
    #              be included in statistic.
    # deleterious: boolean vector of length equal to rows of variants; TRUE if variant is
    #              deleterious, FALSE if it is silent.  Variants with deleterious == NA are
    #              ignored as if loci == FALSE.
    # consider_AA: consider hom alt genotypes in the statistic.
    # consider_RA: consider het genotypes in the statistic.
    # boot:        if TRUE, calculate statistic on a bootstrap sample of the data; if FALSE,
    #              calculate on variants as supplied.

    # Subset to variants of interest
    loci = loci & !is.na(deleterious)
    deleterious = deleterious[loci]
    variants = variants[loci,]

    nRR_test = variants[,paste("nRR_", test_cohort, sep = "")]
    nRA_test = variants[,paste("nRA_", test_cohort, sep = "")]
    nAA_test = variants[,paste("nAA_", test_cohort, sep = "")]
    nmissing_test = variants[,paste("nmissing_", test_cohort, sep = "")]

    nRR_ref = variants[,paste("nRR_", ref_cohort, sep = "")]
    nRA_ref = variants[,paste("nRA_", ref_cohort, sep = "")]
    nAA_ref = variants[,paste("nAA_", ref_cohort, sep = "")]
    nmissing_ref = variants[,paste("nmissing_", ref_cohort, sep = "")]

    # Bootstrap resample, if required
    # Order matters here.  If we resample variants first, then samples next, the
    # 1/3 or so of duplicated variants will likely have different frequencies
    # after the sample drawing.  If we draw samples first, this will not occur.
    # Which is more reasonable?  The point of the variant resampling is to
    # try and capture the uncertainty in the statistic that is due to the arbitrary
    # choice of comparator loci.  Different comparator loci would be unlikely to
    # have identical alt allele frequencies as we are assuming no linkage.
    # The resampling mechanism that best preserves this property is variants then
    # samples.
    if (boot_variants)
    {
        # Resample variants
        resample = sample.int(nrow(variants), replace = TRUE)

        nRR_test = nRR_test[resample]
        nRA_test = nRA_test[resample]
        nAA_test = nAA_test[resample]
        nmissing_test = nmissing_test[resample]
        nRR_ref = nRR_ref[resample]
        nRA_ref = nRA_ref[resample]
        nAA_ref = nAA_ref[resample]
        nmissing_ref = nmissing_ref[resample]

        deleterious = deleterious[resample]
    }

    if (boot_samples)
    {
        # Resample samples
        n_test_redrawn = redraw_samples(nRR_test, nRA_test, nAA_test, nmissing_test)
        nRR_test = n_test_redrawn$RR
        nRA_test = n_test_redrawn$RA
        nAA_test = n_test_redrawn$AA
        nmissing_test = n_test_redrawn$missing

        n_ref_redrawn = redraw_samples(nRR_ref, nRA_ref, nAA_ref, nmissing_ref)
        nRR_ref = n_ref_redrawn$RR
        nRA_ref = n_ref_redrawn$RA
        nAA_ref = n_ref_redrawn$AA
        nmissing_ref = n_ref_redrawn$missing
    }

    # Calculate statistic
    n_test = nAA_test + nRA_test + nRR_test
    f_test_RA = nRA_test / n_test
    f_test_AA = nAA_test / n_test
    n_ref = nAA_ref + nRA_ref + nRR_ref
    f_ref_RA = nRA_ref / n_ref
    f_ref_AA = nAA_ref / n_ref
    test_deleterious = consider_AA*sum(f_test_AA*deleterious, na.rm = TRUE)    + consider_RA*sum(f_test_RA*deleterious, na.rm = TRUE)
    test_silent      = consider_AA*sum(f_test_AA*(!deleterious), na.rm = TRUE) + consider_RA*sum(f_test_RA*(!deleterious), na.rm = TRUE)
    ref_deleterious  = consider_AA*sum(f_ref_AA*deleterious, na.rm = TRUE)     + consider_RA*sum(f_ref_RA*deleterious, na.rm = TRUE)
    ref_silent       = consider_AA*sum(f_ref_AA*(!deleterious), na.rm = TRUE)  + consider_RA*sum(f_ref_RA*(!deleterious), na.rm = TRUE)

    test_rate_deleterious = test_deleterious / (test_deleterious + test_silent)
    ref_rate_deleterious = ref_deleterious / (ref_deleterious + ref_silent)

    test_odds_deleterious = test_rate_deleterious / (1 - test_rate_deleterious)
    ref_odds_deleterious = ref_rate_deleterious / (1 - ref_rate_deleterious)

    # The odds are roughly equivalent to dN/dS ratios, differing in the following aspects:
    # * They are based on average rates across all individuals in a cohort, for the 
    #   genes under test only.
    # * dN and dS are the rate of "expressible" variation, where a variant is only 
    #   considered if it would be predicted to be expressed (eg het genotypes for a 
    #   recessive disorder are not considered).

    # The overall statistic can be interpreted as akin to a difference in log(dN/dS).
    # A value > 0 indicates that the test cohort has a higher dN:dS than the reference.

    stat = log(test_odds_deleterious) - log(ref_odds_deleterious)

    stat
}



enrichment_tests = expand.grid(
    test_cohort = c("mgrb_orig", "mgrb_somfilt", "swegen"),
    ref_cohort = c("swegen", "gnomad"),
    genelist = c("ACMG_ADSD", "ACMG_ADSD_not_somatic"),
    deleterious = c("nonsense", "missense"),
    variants = c("SNV", "all"),
    alpha = 0.05, B = 4000, seed = 314159, stringsAsFactors = FALSE)
enrichment_tests = enrichment_tests[enrichment_tests$test_cohort != enrichment_tests$ref_cohort,]

library(plyr)
library(doParallel)
registerDoParallel(24)


enrichment_tests = ddply(enrichment_tests, colnames(enrichment_tests), function(test_params) {
    set.seed(test_params$seed)

    if (test_params$genelist == "ACMG_ADSD")
    {
        loci = variants$genelist_ACMG_ADSD & variants$rare_all
        consider_AA = TRUE
        consider_RA = TRUE
    }
    else if (test_params$genelist == "ACMG_AR")
    {
        loci = variants$genelist_ACMG_AR & variants$rare_all
        consider_AA = TRUE
        consider_RA = FALSE
    }
    else if (test_params$genelist == "ACMG_ADSD_not_somatic")
    {
        loci = variants$genelist_ACMG_ADSD_not_somatic & variants$rare_all
        consider_AA = TRUE
        consider_RA = TRUE
    }
    else if (test_params$genelist == "ACMG_AR_not_somatic")
    {
        loci = variants$genelist_ACMG_AR_not_somatic & variants$rare_all
        consider_AA = TRUE
        consider_RA = FALSE
    }
    else
        stop("Invalid value for genelist")

    if (test_params$variants == "SNV")
        loci = loci & (variants$variant_class == "SNV")
    else if (test_params$variants == "indel")
        loci = loci & (variants$variant_class %in% c("small_insertion", "small_deletion", "insdel"))
    else if (test_params$variants == "all")
        loci = loci & (variants$variant_class %in% c("SNV", "small_insertion", "small_deletion", "insdel"))
    else
        stop("Invalid value for variant type")

    if (test_params$deleterious == "nonsense")
        deleterious = variants$trunc_vs_syn == "truncating"
    else if (test_params$deleterious == "missense")
        deleterious = variants$nonsyn_vs_syn == "nonsynonymous"
    else
        stop("Invalid value for deleterious")

    stat = enrichment_statistic(variants, loci = loci, deleterious = deleterious, test_cohort = test_params$test_cohort, ref_cohort = test_params$ref_cohort,
        consider_AA = consider_AA, consider_RA = consider_RA, boot_variants = FALSE, boot_samples = FALSE)

    if (is.na(stat))
        boot = NA
    else
        boot = replicate(test_params$B, enrichment_statistic(variants, loci = loci, deleterious = deleterious, test_cohort = test_params$test_cohort, ref_cohort = test_params$ref_cohort,
            consider_AA = consider_AA, consider_RA = consider_RA, boot_variants = TRUE, boot_samples = TRUE))

    # BC interval (http://www.jstor.org/stable/3314608)
    if (any(is.na(boot)) | is.na(stat))
    {
        z0 = NA
        bc_lower = NA
        bc_upper = NA
    }
    else
    {
        z0 = qnorm(mean(boot < stat))
        p_lower = pnorm(2*z0 + qnorm(test_params$alpha/2))
        p_upper = pnorm(2*z0 + qnorm(1-test_params$alpha/2))
        bc_lower = quantile(boot, p_lower)
        bc_upper = quantile(boot, p_upper)
    }

    test_params$statistic = stat
    test_params$z0 = z0
    test_params$boot_lcl = bc_lower
    test_params$boot_median = median(boot)
    test_params$boot_ucl = bc_upper

    test_params
}, .parallel = TRUE)



enrichment_tests


library(openxlsx)
write.xlsx(variants, file = "../MGRB_GnomAD_SweGen_variant_burden.xlsx")


