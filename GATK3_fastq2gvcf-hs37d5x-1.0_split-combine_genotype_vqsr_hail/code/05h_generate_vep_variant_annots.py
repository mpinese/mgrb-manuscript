#!./bin/pyhail.sh
import hail
import subprocess
from hail.expr import *


hc = hail.HailContext(log = 'log/05h_generate_vep_variant_annots.log', tmp_dir = 'tmp/hail')


# Prepare split variants for input to VEP
split_vds = (hc
    .read('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.vds')
    .drop_samples()
    .split_multi()
    .repartition(10000)
    .annotate_variants_expr('va = {rsid:str(v)}')
)
split_vds.write('tmp/05_mgrb_split_variants.vds')

split_vds.export_vcf('tmp/05_mgrb_split_variants.vcf.bgz')


# Call VEP
subprocess.call('bash ./05hh_vep_helper.sh', shell=True)


# Import the VEP results back in to Hail
vep_result = hc.import_table(paths='tmp/05_mgrb_split_variants.vep.tab', comment='##', missing='-', types={
    '#Uploaded_variation': TVariant(),
    'Location': TString(),
    'Allele': TString(),
    'Gene': TString(),
    'Feature': TString(),
    'Feature_type': TString(),
    'Consequence': TString(),
    'cDNA_position': TString(),
    'CDS_position': TString(),
    'Protein_position': TString(),
    'Amino_acids': TString(),
    'Codons': TString(),
    'Existing_variation': TString(),
    'ALLELE_NUM': TInt(),
    'IMPACT': TString(),
    'DISTANCE': TInt(),
    'STRAND': TInt(),
    'FLAGS': TString(),
    'VARIANT_CLASS': TString(),
    'MINIMISED': TInt(),
    'SYMBOL': TString(),
    'SYMBOL_SOURCE': TString(),
    'HGNC_ID': TInt(),
    'BIOTYPE': TString(),
    'CANONICAL': TString(),
    'TSL': TString(),
    'APPRIS': TString(),
    'CCDS': TString(),
    'ENSP': TString(),
    'SWISSPROT': TString(),
    'TREMBL': TString(),
    'UNIPARC': TString(),
    'GENE_PHENO': TInt(),
    'SIFT': TString(),
    'PolyPhen': TString(),
    'EXON': TString(),
    'INTRON': TString(),
    'DOMAINS': TString(),
    'HGVSc': TString(),
    'HGVSp': TString(),
    'HGVS_OFFSET': TInt(),
    'AF': TString(),
    'AFR_AF': TFloat(),
    'AMR_AF': TFloat(),
    'EAS_AF': TFloat(),
    'EUR_AF': TFloat(),
    'SAS_AF': TFloat(),
    'AA_AF': TFloat(),
    'EA_AF': TFloat(),
    'gnomAD_AF': TFloat(),
    'gnomAD_AFR_AF': TFloat(),
    'gnomAD_AMR_AF': TFloat(),
    'gnomAD_ASJ_AF': TFloat(),
    'gnomAD_EAS_AF': TFloat(),
    'gnomAD_FIN_AF': TFloat(),
    'gnomAD_NFE_AF': TFloat(),
    'gnomAD_OTH_AF': TFloat(),
    'gnomAD_SAS_AF': TFloat(),
    'MAX_AF': TFloat(),
    'MAX_AF_POPS': TString(),
    'CLIN_SIG': TString(),
    'SOMATIC': TString(),
    'PHENO': TString(),
    'PUBMED': TString(),
    'MOTIF_NAME': TString(),
    'MOTIF_POS': TInt(),
    'HIGH_INF_POS': TString(),
    'MOTIF_SCORE_CHANGE': TFloat()})

vep_result = (vep_result
    .annotate('''
        Consequences = Consequence.split(","),
        Domains = DOMAINS.split(","),
        ClinVar = CLIN_SIG.split(','),
        PolyPhen_call = PolyPhen.replace("\\\\([0-9.]+\\\\)", ""),
        SIFT_call = SIFT.replace("\\\\([0-9.]+\\\\)", ""),
        PUBMED_list = PUBMED.split(",").map(s => s.toInt()),
        AF_fixed = AF.replace(",.*", "").toFloat(),
        MAX_AF_POPS_list = MAX_AF_POPS.split(",")
    ''')
)

vep_result = (vep_result
    .drop(['Location', 'Allele', 'ALLELE_NUM', 'FLAGS', 'MINIMISED', 'TSL', 'APPRIS', 'GENE_PHENO', 'TREMBL', 'SOMATIC', 'PHENO', 'Consequence', 'DOMAINS', 'CLIN_SIG', 'SIFT', 'PolyPhen'])
    .rename({'#Uploaded_variation': 'Variant'})
    .key_by('Variant')
)


# Annotate the split VDS with the VEP result
split_vds = split_vds.annotate_variants_table(table=vep_result, expr='''
    va.vep = {
        consequences: {
            feature: {
                id: table.Feature,
                type: table.Feature_type,
                biotype: table.BIOTYPE,
                strand: table.STRAND,
                canonical: table.CANONICAL
            },
            position_in_feature: {
                cDNA: table.cDNA_position,
                CDS: table.CDS_position,
                protein: table.Protein_position,
                exon: table.EXON,
                intron: table.INTRON,
                overlapping_domains: table.Domains
            },
            variant_class: table.VARIANT_CLASS,
            distance_from_feature: table.DISTANCE,
            consequences: table.Consequences,
            impact_class: table.IMPACT,
            codon_change: table.Codons,
            residue_change: table.Amino_acids,
            motif_disruption: {
                motif: table.MOTIF_NAME,
                motif_pos: table.MOTIF_POS,
                high_inf_pos: table.HIGH_INF_POS,
                motif_score_change: table.MOTIF_SCORE_CHANGE
            }
        },
        predictions: {
            ClinVar: table.ClinVar,
            PolyPhen: table.PolyPhen_call,
            SIFT: table.SIFT_call
        },
        HGVS: {
           HGVSc: table.HGVSc,
           HGVSp: table.HGVSp,
           HGVS_offset: table.HGVS_OFFSET
        },
        population_frequencies: {
            TGP: {
                AF: table.AF_fixed,
                AFR_AF: table.AFR_AF,
                AMR_AF: table.AMR_AF,
                EAS_AF: table.EAS_AF,
                EUR_AF: table.EUR_AF,
                SAS_AF: table.SAS_AF
            },
            ESP: {
                AA_AF: table.AA_AF,
                EA_AF: table.EA_AF
            },
            GnomAD: {
                AF: table.gnomAD_AF,
                AFR_AF: table.gnomAD_AFR_AF,
                AMR_AF: table.gnomAD_AMR_AF,
                ASJ_AF: table.gnomAD_ASJ_AF,
                EAS_AF: table.gnomAD_EAS_AF,
                FIN_AF: table.gnomAD_FIN_AF,
                NFE_AF: table.gnomAD_NFE_AF,
                OTH_AF: table.gnomAD_OTH_AF,
                SAS_AF: table.gnomAD_SAS_AF
            },
            MAX_AF: table.MAX_AF,
            MAX_AF_POPS: table.MAX_AF_POPS_list
        },
        cross_references: {
            Ensembl_gene: table.Gene,
            Ensembl_protein: table.ENSP,
            PUBMED: table.PUBMED_list,
            CCDS: table.CCDS,
            SWISSPROT: table.SWISSPROT,
            UNIPARC: table.UNIPARC,
            symbol: table.SYMBOL,
            symbol_source: table.SYMBOL_SOURCE,
            HGNC: table.HGNC_ID,
            existing_variation: table.Existing_variation
        }
    }''')


# Write the annotated file for merging with the MGRB multiallelic file.
split_vds.write('tmp/05_mgrb_split_variants.vep.vds')


