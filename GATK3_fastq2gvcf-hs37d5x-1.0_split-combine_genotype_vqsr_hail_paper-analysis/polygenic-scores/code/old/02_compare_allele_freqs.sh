#!/bin/bash

# ../gwas_catalog_v1.0.1-associations_e91_r2018-03-13.tsv downloaded from
# https://www.ebi.ac.uk/gwas/docs/file-downloads, description
# "All associations v1.0.1 - with added ontology annotations and GWAS Catalog study accession numbers"
# on 16/3/2018


mkdir -p tmp
head -n1 ../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.popstar.dosages > tmp/head
LC_ALL=C grep -F -f <(cut -f2 ../manual_polygenic_scores.GnomAD_NFE_AFs.models | tail -n+2) ../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.popstar.dosages > tmp/dosages
cat tmp/head tmp/dosages > ../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.popstar.modelonly.dosages
python3 popstar2tab.py < ../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.popstar.modelonly.dosages > ../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.popstar.modelonly.dosages.tab
