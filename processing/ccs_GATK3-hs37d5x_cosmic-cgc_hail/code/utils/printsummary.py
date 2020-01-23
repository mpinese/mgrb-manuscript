#!./bin/pyhail.sh
import hail
import sys
import pprint

hc = hail.HailContext(tmp_dir = 'tmp/hail')
vds = hc.read(sys.argv[1])

pprint.pprint(vds.summarize())

print('Singletons: {}'.format(vds.filter_variants_expr('gs.filter(g => g.isCalledNonRef()).count() == 1').count_variants()))
