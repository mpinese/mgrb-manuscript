#!/bin/bash

mkdir -p ../source_data
cd ../source_data
rm -f variant_summary.txt.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
date > variant_summary.txt.gz.download_date
cd -
