#!/usr/bin/env python
import sys

state = 0
# States:
# 0  Reading preamble and header
# 1  Reading data

for line in sys.stdin:
    if state == 0:
        # Preamble or header
        sys.stdout.write(line)
        if line.startswith('#CHROM'):
            # Header
            state = 1
    else:
        # Data
        lineparts = line.rstrip().split('\t')
        ref = lineparts[3]
        alts = lineparts[4].split(',')
        data = dict(zip(lineparts[8].split(':'), lineparts[9].split(':')))

        gt = data['GT'].split('/')

        if gt[0] == '.' or gt[0] == '0':
            allele_1 = ref
        else:
            allele_1 = alts[int(gt[0])-1]

        if gt[1] == '.' or gt[1] == '0':
            allele_2 = ref
        else:
            allele_2 = alts[int(gt[1])-1]

        # Skip unwanted sites
        if (allele_1 == ref and allele_2 == ref) or (allele_1 == ref and allele_2 == '*') or (allele_1 == '*' and allele_2 == ref) or (allele_1 == '*' and allele_2 == '*'):
            continue

        # Recode the alts and gt to a minimal representation
        new_alts = list(set([allele for allele in (allele_1, allele_2) if allele != ref]))
        new_gt = sorted([0 if allele == ref else new_alts.index(allele) + 1 for allele in (allele_1, allele_2)])

        lineparts[4] = ','.join(new_alts)
        lineparts[8] = 'GT'
        lineparts[9] = '/'.join([str(gt) for gt in new_gt])

        sys.stdout.write('\t'.join(lineparts) + '\n')
 
