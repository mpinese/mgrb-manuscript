#!./bin/pyhail.sh
import hail
import sys
import pprint

hc = hail.HailContext(tmp_dir = 'tmp/hail')
vds = hc.read(sys.argv[1])

print(vds.num_samples)

