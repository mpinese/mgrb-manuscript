set -euo pipefail

mkdir -p ../source_data

cd ../source_data
wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz
tabix HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz
cd ../code
