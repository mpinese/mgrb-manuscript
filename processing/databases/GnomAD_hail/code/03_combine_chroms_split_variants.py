#!./bin/pyhail.sh
import hail

hc = hail.HailContext(log = 'log/03_combine_chroms_split_variants.log', tmp_dir = 'tmp/hail')

vds_gnomad = hc.read(['../gnomad.genomes.r2.0.1.sites.X.simplified.vds', '../gnomad.genomes.r2.0.1.sites.autosomes.simplified.vds'])

# Code for sex-specific frequency annotations, which are unique to the X file:
# These have been removed in 02_simplify.py to permit loading of the X and 
# autosome files together.
#
#        va.info.AF_AFR_Male = va.info.AF_AFR_Male[va.aIndex - 1],
#        va.info.AF_AMR_Male = va.info.AF_AMR_Male[va.aIndex - 1],
#        va.info.AF_ASJ_Male = va.info.AF_ASJ_Male[va.aIndex - 1],
#        va.info.AF_EAS_Male = va.info.AF_EAS_Male[va.aIndex - 1],
#        va.info.AF_FIN_Male = va.info.AF_FIN_Male[va.aIndex - 1],
#        va.info.AF_NFE_Male = va.info.AF_NFE_Male[va.aIndex - 1],
#        va.info.AF_OTH_Male = va.info.AF_OTH_Male[va.aIndex - 1],
#        va.info.AF_Male = va.info.AF_Male[va.aIndex - 1],
#        va.info.AF_AFR_Female = va.info.AF_AFR_Female[va.aIndex - 1],
#        va.info.AF_AMR_Female = va.info.AF_AMR_Female[va.aIndex - 1],
#        va.info.AF_ASJ_Female = va.info.AF_ASJ_Female[va.aIndex - 1],
#        va.info.AF_EAS_Female = va.info.AF_EAS_Female[va.aIndex - 1],
#        va.info.AF_FIN_Female = va.info.AF_FIN_Female[va.aIndex - 1],
#        va.info.AF_NFE_Female = va.info.AF_NFE_Female[va.aIndex - 1],
#        va.info.AF_OTH_Female = va.info.AF_OTH_Female[va.aIndex - 1],
#        va.info.AF_Female = va.info.AF_Female[va.aIndex - 1],
#        va.info.AC_AFR_Male = va.info.AC_AFR_Male[va.aIndex - 1],
#        va.info.AC_AMR_Male = va.info.AC_AMR_Male[va.aIndex - 1],
#        va.info.AC_ASJ_Male = va.info.AC_ASJ_Male[va.aIndex - 1],
#        va.info.AC_EAS_Male = va.info.AC_EAS_Male[va.aIndex - 1],
#        va.info.AC_FIN_Male = va.info.AC_FIN_Male[va.aIndex - 1],
#        va.info.AC_NFE_Male = va.info.AC_NFE_Male[va.aIndex - 1],
#        va.info.AC_OTH_Male = va.info.AC_OTH_Male[va.aIndex - 1],
#        va.info.AC_Male = va.info.AC_Male[va.aIndex - 1],
#        va.info.AC_AFR_Female = va.info.AC_AFR_Female[va.aIndex - 1],
#        va.info.AC_AMR_Female = va.info.AC_AMR_Female[va.aIndex - 1],
#        va.info.AC_ASJ_Female = va.info.AC_ASJ_Female[va.aIndex - 1],
#        va.info.AC_EAS_Female = va.info.AC_EAS_Female[va.aIndex - 1],
#        va.info.AC_FIN_Female = va.info.AC_FIN_Female[va.aIndex - 1],
#        va.info.AC_NFE_Female = va.info.AC_NFE_Female[va.aIndex - 1],
#        va.info.AC_OTH_Female = va.info.AC_OTH_Female[va.aIndex - 1],
#        va.info.AC_Female = va.info.AC_Female[va.aIndex - 1],
#
#        AF_AFR_Male, AF_AMR_Male, AF_ASJ_Male, AF_EAS_Male, AF_FIN_Male, AF_NFE_Male, AF_OTH_Male, AF_Male, 
#        AF_AFR_Female, AF_AMR_Female, AF_ASJ_Female, AF_EAS_Female, AF_FIN_Female, AF_NFE_Female, AF_OTH_Female, AF_Female,
#        AC_AFR_Male, AC_AMR_Male, AC_ASJ_Male, AC_EAS_Male, AC_FIN_Male, AC_NFE_Male, AC_OTH_Male, AC_Male,
#        AC_AFR_Female, AC_AMR_Female, AC_ASJ_Female, AC_EAS_Female, AC_FIN_Female, AC_NFE_Female, AC_OTH_Female, AC_Female,
#        AN_AFR_Male, AN_AMR_Male, AN_ASJ_Male, AN_EAS_Male, AN_FIN_Male, AN_NFE_Male, AN_OTH_Male, AN_Male,
#        AN_AFR_Female, AN_AMR_Female, AN_ASJ_Female, AN_EAS_Female, AN_FIN_Female, AN_NFE_Female, AN_OTH_Female, AN_Female)


vds_gnomad_split = (vds_gnomad
    .split_multi()
    .annotate_variants_expr('''
        va.info.AF_AFR = if (va.wasSplit) va.info.AF_AFR[va.aIndex - 1] else va.info.AF_AFR[0],
        va.info.AF_AMR = if (va.wasSplit) va.info.AF_AMR[va.aIndex - 1] else va.info.AF_AMR[0],
        va.info.AF_ASJ = if (va.wasSplit) va.info.AF_ASJ[va.aIndex - 1] else va.info.AF_ASJ[0],
        va.info.AF_EAS = if (va.wasSplit) va.info.AF_EAS[va.aIndex - 1] else va.info.AF_EAS[0],
        va.info.AF_FIN = if (va.wasSplit) va.info.AF_FIN[va.aIndex - 1] else va.info.AF_FIN[0],
        va.info.AF_NFE = if (va.wasSplit) va.info.AF_NFE[va.aIndex - 1] else va.info.AF_NFE[0],
        va.info.AF_OTH = if (va.wasSplit) va.info.AF_OTH[va.aIndex - 1] else va.info.AF_OTH[0],
        va.info.AF = if (va.wasSplit) va.info.AF[va.aIndex - 1] else va.info.AF[0],

        va.info.AC_AFR = if (va.wasSplit) va.info.AC_AFR[va.aIndex - 1] else va.info.AC_AFR[0],
        va.info.AC_AMR = if (va.wasSplit) va.info.AC_AMR[va.aIndex - 1] else va.info.AC_AMR[0],
        va.info.AC_ASJ = if (va.wasSplit) va.info.AC_ASJ[va.aIndex - 1] else va.info.AC_ASJ[0],
        va.info.AC_EAS = if (va.wasSplit) va.info.AC_EAS[va.aIndex - 1] else va.info.AC_EAS[0],
        va.info.AC_FIN = if (va.wasSplit) va.info.AC_FIN[va.aIndex - 1] else va.info.AC_FIN[0],
        va.info.AC_NFE = if (va.wasSplit) va.info.AC_NFE[va.aIndex - 1] else va.info.AC_NFE[0],
        va.info.AC_OTH = if (va.wasSplit) va.info.AC_OTH[va.aIndex - 1] else va.info.AC_OTH[0],
        va.info.AC = if (va.wasSplit) va.info.AC[va.aIndex - 1] else va.info.AC[0]''')
    .annotate_variants_expr('''va.gnomad = select(va.info,
        AF_AFR, AF_AMR, AF_ASJ, AF_EAS, AF_FIN, AF_NFE, AF_OTH, AF, 
        AC_AFR, AC_AMR, AC_ASJ, AC_EAS, AC_FIN, AC_NFE, AC_OTH, AC, 
        AN_AFR, AN_AMR, AN_ASJ, AN_EAS, AN_FIN, AN_NFE, AN_OTH, AN)''')
    .annotate_variants_expr('va = select(va, gnomad, wasSplit, aIndex)')
    .naive_coalesce(1000)
    .min_rep()
)
vds_gnomad_split.write('../gnomad.genomes.r2.0.1.sites.combined.split.minrep.vds')
