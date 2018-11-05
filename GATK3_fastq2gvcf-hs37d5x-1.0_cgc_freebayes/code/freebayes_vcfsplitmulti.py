#!/usr/bin/env python3

# Output: VCF with multi-allelics split


def process_vcf(infile, outfile):
    split_i = 1
    for line in infile:
        if line.startswith('#CHROM'):
            outfile.write('##INFO=<ID=SPLIT,Number=1,Type=Integer,Description="Index of the split group for a split allele, or 0 if not split.">\n')
        if line.startswith('#'):
            outfile.write(line)
            continue
        lineparts = line.rstrip().split('\t')
        chrom, pos, vid, ref, alt_str, qual, filt, info_str, fmt_str, data_str = lineparts

        data = dict(zip(fmt_str.split(':'), data_str.split(':')))
        info = dict([x.split('=') for x in info_str.split(';')])

        alts = alt_str.split(',')
        aos = data['AO'].split(',')
        qas = data['QA'].split(',')
        rpls = info['RPL'].split(',')
        rprs = info['RPR'].split(',')
        rpps = info['RPP'].split(',')
        safs = info['SAF'].split(',')
        sars = info['SAR'].split(',')
        saps = info['SAP'].split(',')
        epps = info['EPP'].split(',')
        paireds = info['PAIRED'].split(',')
        mqms = info['MQM'].split(',')
        ro = int(data['RO'])
        dp = int(data['DP'])

        if len(alts) > 1:
            info['SPLIT'] = str(split_i)
            split_i += 1
        else:
            info['SPLIT'] = '0'

        for alt, ao, qa, rpl, rpr, rpp, saf, sar, sap, epp, paired, mqm in zip(alts, aos, qas, rpls, rprs, rpps, safs, sars, saps, epps, paireds, mqms):
            info_i = info.copy()
            info_i['RPL'] = rpl
            info_i['RPR'] = rpr
            info_i['RPP'] = rpp
            info_i['SAF'] = saf
            info_i['SAR'] = sar
            info_i['SAP'] = sap
            info_i['EPP'] = epp
            info_i['PAIRED'] = paired
            info_i['MQM'] = mqm
            outfile.write('\t'.join([chrom, pos, vid, ref, alt, qual, filt, ';'.join(['='.join(x) for x in info_i.items()]), 'DP:RO:QR:AO:QA', ':'.join([data['DP'], data['RO'], data['QR'], ao, qa])]) + '\n')


if __name__ == '__main__':
    import sys
    process_vcf(sys.stdin, sys.stdout)

