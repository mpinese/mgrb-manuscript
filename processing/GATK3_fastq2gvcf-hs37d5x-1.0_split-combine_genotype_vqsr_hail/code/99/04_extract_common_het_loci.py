#!../../software/pyhail.sh
import hail
from hail.representation import Interval

hc = hail.HailContext(log = 'log/02_export_split_variants.log', tmp_dir = 'tmp/hail')

vds = hc.read('../MGRB.phase2.tier12.match.vqsr.minrep.vds')

# Stringently filter variants down to common autosomal SNPs with no evidence of HWE deviation.
vds = (vds
    .filter_intervals([Interval.parse('X'), Interval.parse('Y'), Interval.parse('MT')], keep=False)
    .filter_multi()
    .variant_qc()
    .filter_variants_expr('''
        v.altAllele.isSNP &&
        va.qc.callRate >= 0.99 && 
        va.qc.pHWE > 1e-3 && 
        va.qc.dpMean >= 20 && va.qc.dpMean <= 60 && 
        va.qc.dpStDev < 8 && 
        va.qc.rHeterozygosity >= 0.2 && va.qc.rHeterozygosity <= 0.8''')
)

# Drop samples with poor metrics on these filtered variants.
vds = (vds
    .sample_qc()
    .filter_samples_expr('sa.qc.callRate >= 0.985')
)

vds = vds.repartition(28*20, shuffle=False)

#vds.write('../MGRB.phase2.tier12.match.vqsr.minrep.commonhets.vds')

vds.export_genotypes('../MGRB.phase2.tier12.match.vqsr.minrep.commonhets.tsv', 'SAMPLE=s, CHROM=v.contig, POS=v.start, DP=g.dp, AD=g.ad[1]', export_ref=False, export_missing=False)


outfiles = {}
infile = open('../MGRB.phase2.tier12.match.vqsr.minrep.commonhets.tsv', 'rt')
next(infile)
for line in infile:
    sample, chrom, pos, dp, ad = line.rstrip().split('\t')

    if dp == ad or ad == '0':
        continue

    if sample not in outfiles:
        outfiles[sample] = open('../MGRB.phase2.tier12.match.vqsr.minrep.commonhets.{}.tsv'.format(sample), 'wt')

    outfiles[sample].write('\t'.join((chrom, pos, dp, ad)) + '\n')
