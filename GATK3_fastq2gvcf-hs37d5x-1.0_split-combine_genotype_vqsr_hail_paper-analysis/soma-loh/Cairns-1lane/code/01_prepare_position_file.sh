#!/bin/bash
set -euo pipefail

sed 's/:/ /g' MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.commonhets.variant_list | cut -d' ' -f1-2 | sort -k1,1 -k2,2n > 01_positions.list
