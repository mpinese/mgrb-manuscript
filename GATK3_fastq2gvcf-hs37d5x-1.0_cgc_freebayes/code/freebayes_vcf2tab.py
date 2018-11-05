#!/usr/bin/env python3
import warnings


# Output: TSV with columns:
# sample        
# chrom         
# pos           
# ref           
# alt           
# dp            Read Depth
# rd            Number of observations for the ref allele
# ad            Number of observations for the given alt allele
# ado           Number of observations of any other alt alleles (if multiallelic)
# vaf           VAF of the alt allele relative to all reads (ad/dp)
# vafo          VAF of the alt allele discounting any other alt alleles (ad/(ad+rd))
# qr            Sum of quality of the reference observations
# qa            Sum of quality of the alternate observations
# split         Index of the split group if this locus is a split multiallelic, else 0
# mqm           Mean mapping quality of observed alternate alleles
# mqmr          Mean mapping quality of observed reference alleles
# paired        Proportion of observed alternate alleles which are supported by properly paired read fragments
# pairedr       Proportion of observed reference alleles which are supported by properly paired read fragments
# epp           End Placement Probability: Phred-scaled upper-bounds estimate of the probability of observing the deviation between EL and ER given E(EL/ER) ~ 0.5
# eppr          End Placement Probability for reference observations: Phred-scaled upper-bounds estimate of the probability of observing the deviation between EL and ER given E(EL/ER) ~ 0.5
# saf           Number of alternate observations on the forward strand
# sar           Number of alternate observations on the reverse strand
# sap           Strand balance probability for the alternate allele: Phred-scaled upper-bounds estimate of the probability of observing the deviation between SAF and SAR given E(SAF/SAR) ~ 0.5
# srf           Number of reference observations on the forward strand
# srr           Number of reference observations on the reverse strand
# srp           Strand balance probability for the reference allele: Phred-scaled upper-bounds estimate of the probability of observing the deviation between SRF and SRR given E(SRF/SRR) ~ 0.5
# rpl           Reads Placed Left: number of reads supporting the alternate balanced to the left (5') of the alternate allele
# rpr           Reads Placed Right: number of reads supporting the alternate balanced to the right (3') of the alternate allele
# rpp           Read Placement Probability: Phred-scaled upper-bounds estimate of the probability of observing the deviation between RPL and RPR given E(RPL/RPR) ~ 0.5
# rppr          Read Placement Probability for reference observations: Phred-scaled upper-bounds estimate of the probability of observing the deviation between RPL and RPR given E(RPL/RPR) ~ 0.5


def process_vcf(infile, outfile):
    outfile.write('\t'.join([
        'sample', 'chrom', 'pos', 'ref', 'alt', 'dp', 'rd', 'ad', 'ado', 'vaf', 'vafo', 'qr', 'qa',
        'split', 'mqm', 'mqmr', 'paired', 'pairedr', 'epp', 'eppr', 
        'saf', 'sar', 'sap', 'srf', 'srr', 'srp', 'rpl', 'rpr', 'rpp', 'rppr']) + '\n')
    for line in infile:
        if line.startswith('##'):
            continue
        lineparts = line.rstrip().split('\t')
        if line.startswith('#'):
            assert len(lineparts) == 10
            sampleid = lineparts[9]
            continue
        chrom, pos, vid, ref, alt, qual, filt, info_str, fmt_str, data_str = lineparts

        info = {}
        for field in info_str.split(';'):
            field_parts = field.split('=', 1)
            if len(field_parts) == 1:
                info[field_parts[0]] = True
            else:
                info[field_parts[0]] = field_parts[1]

        data = dict(zip(fmt_str.split(':'), data_str.split(':')))

        if ',' in alt:
            warnings.warn('Skipping multi-allelic locus; split with freebayes_vcfsplitmulti.py first')
            continue

        ao = int(data['AO'])
        qa = data['QA']
        ro = int(data['RO'])
        dp = int(data['DP'])

        oad = dp - ro - ao
        vaf = ao*1.0 / dp
        vafs = ao*1.0 / (ro + ao)
        outfile.write('\t'.join([
            sampleid, chrom, pos, ref, alt, data['DP'], data['RO'], str(ao), str(oad), str(vaf), str(vafs), data['QR'], qa,
            info['SPLIT'], info['MQM'], info['MQMR'], info['PAIRED'], info['PAIREDR'], info['EPP'], info['EPPR'],
            info['SAF'], info['SAR'], info['SAP'], info['SRF'], info['SRR'], info['SRP'], info['RPL'], info['RPR'], info['RPP'], info['RPPR']]) + '\n')


if __name__ == '__main__':
    import sys
    process_vcf(sys.stdin, sys.stdout)

