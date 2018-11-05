#!./bin/pyhail.sh
import hail

hc = hail.HailContext(log = 'log/04_import_to_hail.log', tmp_dir = 'tmp/hail')


# For some reason the direct approach consistently crashes
# with a Spark worker timeout, therefore this ridiculous 
# workaround.  I have no idea why this works.
for chrom in [str(n) for n in range(1, 23)] + ['X', 'Y', 'MT']:
    (hc
        .import_vcf('../normed_data/All_20170710.normed.chrom{}.vcf.bgz'.format(chrom))
        .annotate_variants_expr('va.info = {rsid: va.rsid}')
        .min_rep()
        .repartition(28)
        .write('../dbSNP_150.{}.minrep.vds'.format(chrom))
    )


vds = hc.read(['../dbSNP_150.{}.minrep.vds'.format(chrom) for chrom in [str(n) for n in range(1, 6)]])
vds.repartition(280).write('../dbSNP_150.chrom1-5.minrep.vds')

vds = hc.read(['../dbSNP_150.{}.minrep.vds'.format(chrom) for chrom in [str(n) for n in range(6, 11)]])
vds.repartition(280).write('../dbSNP_150.chrom6-10.minrep.vds')

vds = hc.read(['../dbSNP_150.{}.minrep.vds'.format(chrom) for chrom in [str(n) for n in range(11, 16)]])
vds.repartition(280).write('../dbSNP_150.chrom11-15.minrep.vds')

vds = hc.read(['../dbSNP_150.{}.minrep.vds'.format(chrom) for chrom in [str(n) for n in range(16, 21)]])
vds.repartition(280).write('../dbSNP_150.chrom16-20.minrep.vds')

vds = hc.read(['../dbSNP_150.{}.minrep.vds'.format(chrom) for chrom in [str(n) for n in range(21, 23)]])
vds.repartition(280).write('../dbSNP_150.chrom21-22.minrep.vds')

vds = hc.read(['../dbSNP_150.{}.minrep.vds'.format(chrom) for chrom in ['X', 'Y', 'MT']])
vds.repartition(280).write('../dbSNP_150.chromXYMT.minrep.vds')


hc.read(['../dbSNP_150.chrom1-5.minrep.vds', '../dbSNP_150.chrom6-10.minrep.vds']).write('../dbSNP_150.chrom1-10.minrep.vds')
hc.read(['../dbSNP_150.chrom11-15.minrep.vds', '../dbSNP_150.chrom16-20.minrep.vds']).write('../dbSNP_150.chrom11-20.minrep.vds')


hc.read(['../dbSNP_150.chrom1-10.minrep.vds', '../dbSNP_150.chrom11-20.minrep.vds', '../dbSNP_150.chrom21-22.minrep.vds']).write('../dbSNP_150.autosomes.minrep.vds')

hc.read(['../dbSNP_150.autosomes.minrep.vds', '../dbSNP_150.chromXYMT.minrep.vds']).split_multi().min_rep().repartition(560).write('../dbSNP_150.split.minrep.vds')
