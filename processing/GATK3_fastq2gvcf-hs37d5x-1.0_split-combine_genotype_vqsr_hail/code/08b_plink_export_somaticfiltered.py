#!./bin/pyhail-0.1-latest.sh
import pyspark
import hail
from hail import KeyTable
from hail.representation import Interval

hc = hail.HailContext(log = 'log/08b_plink_export_somaticfiltered.log', tmp_dir = 'tmp/hail')

#vds = hc.read('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.vds')
vds = hc.read('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.vds')

tier1_bed = KeyTable.import_bed('../../locus-annotations/source_data/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed')

# Extract good probably not somatic markers for rare variant comparisons.
# Definition of 'good markers':
#   * In autosomes
#   * In tier 1 regions.
# 
# Definition of probably not somatic:
#   DP > 10 AND
#   (
#     GT != het OR
#     (
#       binomTest(ad, dp, 0.5, "two.sided") >= alpha
#     )
#   )
#
# With alpha = 0.1 we expect a little over 10% of het sites will be falsely lost.
# For alpha = 0.1, the DP > 10 restriction corresponds to effective bounds on VAF
# of [0.2, 0.8].

for chrom in range(1, 23):
    print 'Chromosome %d...' % chrom
    (vds
        .filter_intervals(Interval.parse('%d' % chrom))
        .annotate_variants_table(tier1_bed, root='va.tier1bed')
        .filter_variants_expr('va.tier1bed == true && v.isAutosomal()', keep=True)
        .filter_samples_expr('"^Z" ~ s', keep=False)
        .split_multi()
        .filter_genotypes('''
            g.dp > 10 && 
            (
              (!g.isHet()) || binomTest(g.ad[1], g.dp, 0.5, "two.sided") >= 0.1
            )''', keep=True)
        .min_rep()
        .export_plink('tmp/08b_genotypes_%d' % chrom)
    )
