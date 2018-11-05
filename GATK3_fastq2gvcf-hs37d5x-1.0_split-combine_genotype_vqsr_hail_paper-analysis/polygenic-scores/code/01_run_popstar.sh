#!/bin/bash
set -euo pipefail

./popstar/popstar \
    --dosages=../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.popstar.modelonly.dosages \
    --models=../manual_polygenic_scores.GnomAD_NFE_AFs.models \
    --out=../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.manual_polygenic_scores.GnomAD_NFE_AFs.popstar.externalref.scores \
    --format=complete --iter=0 --ref=external

./popstar/popstar \
    --dosages=../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.popstar.dosages \
    --models=../manual_polygenic_scores.GnomAD_NFE_AFs.models \
    --out=../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.manual_polygenic_scores.GnomAD_NFE_AFs.popstar.externalref.permsummary \
    --format=summary --iter=40000 --ref=external

./popstar/popstar \
    --dosages=../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.popstar.dosages \
    --models=../manual_polygenic_scores.GnomAD_NFE_AFs.models \
    --out=../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.manual_polygenic_scores.GnomAD_NFE_AFs.popstar.internalref.permsummary \
    --format=summary --iter=40000 --ref=internal

./popstar/popstar \
    --dosages=../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.popstar.dosages \
    --models=../manual_polygenic_scores.GnomAD_NFE_AFs.models \
    --out=../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.manual_polygenic_scores.GnomAD_NFE_AFs.popstar.externalref.permscores \
    --format=complete --iter=1000 --ref=external

./popstar/popstar \
    --dosages=../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.popstar.dosages \
    --models=../manual_polygenic_scores.GnomAD_NFE_AFs.models \
    --out=../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.manual_polygenic_scores.GnomAD_NFE_AFs.popstar.internalref.permscores \
    --format=complete --iter=1000 --ref=internal

Rscript 01h_tables_to_rds.R

rm ../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.manual_polygenic_scores.GnomAD_NFE_AFs.popstar.externalref.scores
rm ../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.manual_polygenic_scores.GnomAD_NFE_AFs.popstar.externalref.permsummary
rm ../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.manual_polygenic_scores.GnomAD_NFE_AFs.popstar.internalref.permsummary
rm ../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.manual_polygenic_scores.GnomAD_NFE_AFs.popstar.externalref.permscores
rm ../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.manual_polygenic_scores.GnomAD_NFE_AFs.popstar.internalref.permscores
