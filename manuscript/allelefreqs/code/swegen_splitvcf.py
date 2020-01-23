#!/usr/bin/env python3

# **** WARNING: THIS CODE DROPS MULTI-ALLELIC VARIANTS ****
# **** WARNING: THIS CODE DROPS MULTI-ALLELIC VARIANTS ****
# **** WARNING: THIS CODE DROPS MULTI-ALLELIC VARIANTS ****
# See below

# NOTE: hard-coded for 1000 samples, correct values for autosomes only.
#
# **** WARNING: THIS CODE DROPS MULTI-ALLELIC VARIANTS ****
# This is unavoidable; the data supplied in the SweGen AF VCF are 
# not sufficient to determine genotype frequencies at multi-allelic
# positions.
#
# Due to the dropping of multi-allelics the following old note
# is obsolete, but kept for reference: 
# NOTE: the correct interpretation of star alleles is tricky here.
# For the purposes of this rare variant burden comparison, star alleles
# are downcoded to reference.  This is done for two reasons:
#  * It avoids double-counting variants (which could happen if the star
#    alleles were labelled as anything other than reference, considering
#    that the 'actual' allele leading to the star will also be present).
#  * It avoids distortion of AFs -- consider the case of a very common
#    but inconsequential deletion, and a highly deleterious variant in
#    some individuals in the deleted region.  This deleterious variant
#    has in truth a low population frequency, but a high frequency if the
#    star alleles are not included in calculations.


if __name__ == '__main__':
    import sys

    sys.stdout.write('''##fileformat=VCFv4.2
##INFO=<ID=nRR,Number=1,Type=Integer>
##INFO=<ID=nRA,Number=1,Type=Integer>
##INFO=<ID=nAA,Number=1,Type=Integer>
##INFO=<ID=nmissing,Number=1,Type=Integer>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO''')

    for line in sys.stdin:
        if not line.startswith('#'):
            chrom, pos, vid, ref, alt_string, qual, filt, info_string = line.rstrip().split('\t')

            alt_list = alt_string.split(',')
            if len(alt_list) > 1:
                continue               ### MULTI-ALLELICS DROPPED HERE

            info_dict = {}
            for item in info_string.split(';'):
                if '=' in item:
                    key, value = item.split('=', maxsplit=1)
                else:
                    key = item
                    value = None
                info_dict[key] = value

            nmissing = int(1000 - int(info_dict['AN'])/2)
            nRA_list = [int(x) for x in info_dict['AC_Het'].split(',')]
            nAA_list = [int(int(x)/2) for x in info_dict['AC_Hom'].split(',')]
            nRR = 1000 - nmissing - sum(nAA_list) - sum(nRA_list)

            # Transform star alleles to reference
            if '*' in alt_list:
                star_idx = alt_list.index('*')
                nRR += nAA_list[star_idx] + nRA_list[star_idx]

            # Implicitly split alleles as per Hail's default method, 
            # coding a 1/2 variant as two 0/1 variants.
            for alt, nRA, nAA  in zip(alt_list, nRA_list, nAA_list):
                if alt == '*':
                    # Star alleles have already been transformed to reference
                    # in the above.
                    continue
                sys.stdout.write('\n{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t{qual}\t{filt}\tnRR={nRR};nRA={nRA};nAA={nAA};nmissing={nmissing}'.format(
                    chrom=chrom, pos=pos, vid=vid, ref=ref, alt=alt, qual=qual, filt=filt, nRR=nRR, nRA=nRA, nAA=nAA, nmissing=nmissing))
