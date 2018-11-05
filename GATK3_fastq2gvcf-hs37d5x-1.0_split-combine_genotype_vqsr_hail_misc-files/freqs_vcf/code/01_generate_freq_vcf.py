#!./bin/pyhail.sh
import hail
import os
import os.path

if not os.path.exists('tmp'):
    os.mkdir('tmp')

hc = hail.HailContext(log = 'log/01_generate_freq_vcf.log', tmp_dir = 'tmp/hail')

vds_mgrb = hc.read('../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.vds')

# Split variants
vds_mgrb_split = vds_mgrb.split_multi().min_rep()

(vds_mgrb_split
    .annotate_variants_expr('va = {}')
    .variant_qc()
    .drop_samples()
    .annotate_variants_expr('''va.info = {
        AC:       va.qc.AC,
        AF:       va.qc.AF,
        NS:       va.qc.nCalled,
        callRate: va.qc.callRate,
        nHomRef:  va.qc.nHomRef,
        nHet:     va.qc.nHet,
        nHomVar:  va.qc.nHomVar 
    }''')
    .set_va_attributes('va.info.AC', {'Description': 'Variant allele count', 'Number': '1'})
    .set_va_attributes('va.info.AF', {'Description': 'Variant allele frequency', 'Number': '1'})
    .set_va_attributes('va.info.NS', {'Description': 'Samples with genotype', 'Number': '1'})
    .set_va_attributes('va.info.callRate', {'Description': 'Call rate', 'Number': '1'})
    .set_va_attributes('va.info.nHomRef', {'Description': 'Count of homozygous reference allele individuals', 'Number': '1'})
    .set_va_attributes('va.info.nHet', {'Description': 'Count of heterozygous individuals', 'Number': '1'})
    .set_va_attributes('va.info.nHomVar', {'Description': 'Count of homozygous alternate allele individuals', 'Number': '1'})
    .export_vcf('./tmp/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.freqs.vcf.bgz')
)

