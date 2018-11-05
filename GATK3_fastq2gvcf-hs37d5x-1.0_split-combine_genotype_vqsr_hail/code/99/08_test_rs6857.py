#!../../software/pyhail.sh
import hail

from hail.representation import Interval

hc = hail.HailContext(log = 'log/08_test_rs6857.log', tmp_dir = 'tmp/hail')

vds = hc.read('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.vds')

(vds
    .filter_intervals(Interval.parse('19:45392250-45392260'), keep=True)
    .export_variants('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.rs6857.tsv', 'variant=v, va.*')
)
#    .export_genotypes('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.rs6857.tsv', 'sample=s, variant=v, gt=""+g.gtj()+"/"+g.gtk(), dp=g.dp, dpj=g.ad[g.gtj()], dpk=g.ad[g.gtk()], va.*')
