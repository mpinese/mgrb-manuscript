#!./bin/pyhail.sh
import hail
from hail.expr import *


hc = hail.HailContext(log = 'log/02_import_to_hail.log', tmp_dir = 'tmp/hail')

gavin_table = hc.import_table(
    paths='../source_data/GAVIN_ruleguide_r0.3.tsv',
    comment='##',
    no_header=False,
    delimiter='\t',
    missing='n/a',
    types={
        'Gene': TString(),
        'PathogenicIfCADDScoreGreaterThan': TFloat(),
        'BenignIfCADDScoreLessThan': TFloat(),
        'BenignIfMAFGreaterThan': TFloat(),
        'PathogenicIfImpactEquals': TString(),
        'CalibrationCategory': TString()
    }).key_by('Gene')

gavin_table.write('../GAVIN_ruleguide_r0.3.kt')

'''
################################
## GAVIN Applied Rule Guide r0.3
################################
## 
## Variant can be interpreted by following the columns from left to right.
## This classification scheme was implemented in Java, and used to
## benchmark GAVIN in the paper (Van der Velde et al., using r0.2), see:
## https://github.com/molgenis/molgenis and https://github.com/molgenis/gavin
## 
## Genome-wide rules are used if the gene-specific rules fail to classify.
## These rules are applied as follows:
## 1) If impact equals MODIFIER -> benign,
## 2) if MAF greater than 0.003456145 -> benign,
## 3) if CADD greater than 15 -> pathogenic,
## 4) if CADD less than 15 -> benign.
## 
## Explanation of the gene calibration categories:
## C1 = CADD scores highly significantly predictive for pathogenicity (pval < 0.01).
## C2 = CADD scores significantly predictive for pathogenicity (pval < 0.05).
## C3 = CADD scores may be predictive for pathogenicity (pval > 0.05 but with few samples).
## C4 = CADD scores less predictive for pathogenicity (pval > 0.05 with enough samples).
## C5 = CADD scores less predictive for pathogenicity (population CADD > pathogenic CADD).
## I1 = HIGH impact unique for, thus predictive, for pathogenic variants.
## I2 = MODERATE or HIGH impact unique, thus predictive, for pathogenic variants.
## I3 = LOW or MODERATE or HIGH impact unique, thus predictive, for pathogenic variants.
## T1 = Too few exac variants after filtering with pathogenic 95th percentile MAF.
## T2 = Too few exac variants after filtering with impact distribution.
## N1 = Too few ClinVar variants for calibration at this time.
## N2 = Too few ExAC variants found for calibration.
## 
## For C1 and C2, CADD score thresholds are means of stratified benign and pathogenic variants.
## For C3, C4 and C5, CADD score thresholds are 95th sensitivity/specificity percentiles of stratified benign and pathogenic variants.
## 
Gene	PathogenicIfCADDScoreGreaterThan	BenignIfCADDScoreLessThan	BenignIfMAFGreaterThan	PathogenicIfImpactEquals	CalibrationCategory
NUP107	31.700000000000003	15.27	8.236E-5	n/a	C4
'''
