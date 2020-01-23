#!./bin/pyhail.sh
import hail

hc = hail.HailContext(tmp_dir = 'tmp/hail')
hc.read('../dbSNP_150.minrep.vds').split_multi().min_rep().repartition(560).write('../dbSNP_150.split.minrep.vds')

