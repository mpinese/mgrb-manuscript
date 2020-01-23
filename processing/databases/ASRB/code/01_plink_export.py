#!./bin/pyhail.sh
import pyspark
import hail
from hail import KeyTable
from hail.representation import Interval

hc = hail.HailContext(log = 'log/01_plink_export.log', tmp_dir = 'tmp/hail')

vds = hc.read('../ASRB.WGStier12.split.minrep.vds')

tier1_bed = KeyTable.import_bed('../../../locus-annotations/source_data/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed')

# Extract good markers for rare variant comparisons.
# Definition of 'good markers':
#   * In autosomes
#   * In tier 1 regions.
(vds
    .annotate_variants_table(tier1_bed, root='va.tier1bed')
    .filter_variants_expr('va.tier1bed == true && v.isAutosomal()', keep=True)
    .filter_samples_expr('"^Z" ~ s', keep=False)
    .split_multi()
    .min_rep()
    .export_plink('tmp/01_genotypes')
)
