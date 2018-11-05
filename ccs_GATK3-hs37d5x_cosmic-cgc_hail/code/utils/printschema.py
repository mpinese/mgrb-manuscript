#!./bin/pyhail.sh
import hail
import sys
import pprint

hc = hail.HailContext(tmp_dir = 'tmp/hail')
vds = hc.read(sys.argv[1])

print('Sample schema:')
pprint.pprint(vds.sample_schema)
print('\nVariant schema:')
pprint.pprint(vds.variant_schema)
print('\nGlobal schema:')
pprint.pprint(vds.global_schema)
