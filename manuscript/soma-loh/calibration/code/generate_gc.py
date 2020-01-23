#!/usr/bin/env python
import pysam


def calcgc(seq):
    n_gc = 0
    n_tot = 0
    for base in seq.upper():
        if base == 'G' or base == 'C':
            n_gc += 1
            n_tot += 1
        elif base == 'A' or base == 'T':
            n_tot += 1
    return float(n_gc)/n_tot


if __name__ == '__main__':
    ref = pysam.FastaFile('../../../../../resources/hs37d5x/reference_genome/hs37d5x.fa')

    variant_file = open('../unfiltered/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.commonhets.variant_list', 'rt')

    print('chrom\tpos\tgc100\tgc200\tgc400\tgc600\tgc800')

for line in variant_file:
    chrom, pos, _, _ = line.rstrip().split(':')

    posi = int(pos)
    seq800 = ref.fetch(chrom, start=posi-400, end=posi+400)
    seq600 = seq800[100:-100]
    seq400 = seq600[100:-100]
    seq200 = seq400[100:-100]
    seq100 = seq200[50:-50]

    print chrom + '\t' + pos + '\t' + str(calcgc(seq100)) + '\t' + str(calcgc(seq200)) + '\t' + str(calcgc(seq400)) + '\t' + str(calcgc(seq600)) + '\t' + str(calcgc(seq800))
