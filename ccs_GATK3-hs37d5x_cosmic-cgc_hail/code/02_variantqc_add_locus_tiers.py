#!./bin/pyhail.sh
import hail
from hail.expr import TInt, TVariant
import subprocess
import os

hc = hail.HailContext(log = 'log/02_variantqc_add_locus_tiers.log', tmp_dir = 'tmp/hail')

vds = hc.read('../ccs.cosmic_cgc.minrep.vds')

vds.export_variants('../variants.bed', 'v.contig, v.start-1, v.start+v.ref.length()-1, v')

subprocess.call('sort -T tmp/ -k1,1 -k2,2n ../variants.bed | bgzip > ../variants.sorted.bed.gz', shell=True)
subprocess.call('tabix ../variants.sorted.bed.gz', shell=True)

subprocess.call('gzip -dc ../../locus-annotations/hs37d5x_data/genome.bed.gz | cut -f 1,3 > ../genome.txt', shell=True)

regions = [('../../locus-annotations/regions/bad.depth.bed.gz', '../variants.distance.bad.depth.txt', 'va.locus.badCoverage'),
           ('../../locus-annotations/regions/bad.mdust.bed.gz', '../variants.distance.bad.mdust.txt', 'va.locus.badComplexity'),
           ('../../locus-annotations/regions/bad.rmsk.bed.gz', '../variants.distance.bad.rmsk.txt', 'va.locus.badRepeat'),
           ('../../locus-annotations/regions/bad.wgEncodeCrgMapabilityAlign100mer.bed.gz', '../variants.distance.bad.mappability100.txt', 'va.locus.badMappability'),
           ('../../locus-annotations/regions/bad.wgEncodeExcludable.bed.gz', '../variants.distance.bad.encodeExcludable.txt', 'va.locus.badEncodeExcluded'),
           ('../../locus-annotations/regions/good.HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed.gz', '../variants.distance.good.GiaB_HG001_highConf.txt', 'va.locus.goodGiaBHighConf')]

vds = vds.annotate_variants_expr('va.locus = {}')

for inpath, outpath, root in regions:
    cmd = './bin/bedtools closest -d -g ../genome.txt -a ../variants.sorted.bed.gz -b {inpath} > {outpath}'.format(inpath=inpath, outpath=outpath)
    print(cmd)
    subprocess.call(cmd, shell=True)
    kt = hc.import_table(outpath, no_header=True, types={'f3': TVariant(), 'f7': TInt()}).annotate('v = f3, dist = f7').key_by('v')
    vds = vds.annotate_variants_table(kt, expr='{}_dist = table.dist'.format(root)).annotate_variants_expr('{} = {}_dist == 0'.format(root, root))
#    os.remove(outpath)
	
vds = (vds
    .annotate_variants_expr('va.locus.badPAR = v.inXPar() || v.inYPar()')
    .annotate_variants_expr('''va.locus.tier = 
        if (!("^([0-9]+|X|Y)$" ~ v.contig)) 3
        else (if (va.locus.badCoverage_dist <= 5 || va.locus.badComplexity_dist <= 5 || va.locus.badMappability_dist <= 5 || va.locus.badEncodeExcluded_dist <= 5 || va.locus.badPAR) 3
        else (if (va.locus.goodGiaBHighConf_dist > 0) 2
        else (if (v.contig ~ "^(X|Y)$") 2
        else 1)))''')
)

vds.write('../ccs.cosmic_cgc.minrep.locusannot.vds')
