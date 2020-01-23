#!./bin/pyhail.sh
import hail

hc = hail.HailContext(tmp_dir = 'tmp/hail')

# Merge the MGRB high quality, 1000 genomes, and HRC datasets to get a common 
# set of variants.
vds_mgrb = hc.read('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.vds')

vds_mgrb.export_samples('samples.txt', 's')
