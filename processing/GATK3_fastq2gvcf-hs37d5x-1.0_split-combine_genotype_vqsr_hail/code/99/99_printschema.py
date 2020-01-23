#!../../software/pyhail.sh
import hail
from hail.representation import Interval

hc = hail.HailContext(log = 'log/02_export_split_variants.log', tmp_dir = 'tmp/hail')

vds = hc.read('../MGRB.phase2.tier12.match.vqsr.minrep.vds')

print(vds.sample_schema)

