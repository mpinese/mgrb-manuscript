import re


def parse_header(instream):
    info_formats = {}

    for line in instream:
        lineparts = line.rstrip().split('\t')
        if line.startswith('#CHROM'):
            samples = lineparts[9:]
            break
        elif line.startswith('##INFO='):
            match = re.search('##INFO=<ID=(?P<id>[^,]+),Number=(?P<number>[0-9]+|[ARG.]),Type=(?P<type>Integer|Float|Flag|Character|String),Description="(?P<desc>[^"]*)"', line)
            assert match != None, line
            groups = match.groupdict()
            if groups['type'] == "String" and 'Format' in groups['desc']:
                fmt_match = re.search('Format: +([^=]+=)?(?P<format>[^",]+)', groups['desc'])
                assert fmt_match != None, groups['desc']
                info_formats[groups['id']] = fmt_match.group('format').split('|')

    return samples, info_formats


def missing_to_none(value):
    if value == '.':
        return None
    else:
        return value


def process_data(instream, samples, info_formats, outstream):
#    outstream.write('Sample\tChrom\tPos\tID\tRef\tGenotype\tAlt\tAltDepth\tTotalDepth\tMGRB_AC\tMGRB_AF')
    outstream.write('Sample\tChrom\tPos\tID\tRef\tGenotype\tAlt\tAltDepth\tTotalDepth\tTier')
    for vep_field in info_formats['CSQ']:
        outstream.write('\tVEP.{}'.format(vep_field))
#    for eigen_field in info_formats['Eigen']:
#        outstream.write('\tEigen.{}'.format(eigen_field))
#    for cato_field in info_formats['Cato3']:
#        outstream.write('\tCato3.{}'.format(cato_field))
#    for iarc_field in info_formats['IARCTP53']:
#        outstream.write('\tIARCTP53.{}'.format(iarc_field))
    outstream.write('\n')

    for line in instream:
        lineparts = line.rstrip().split('\t')

        chrom = lineparts[0]
        pos = lineparts[1]
        vid = missing_to_none(lineparts[2])
        ref = lineparts[3]
        alts = [missing_to_none(a) for a in lineparts[4].split(',')]
        alleles = [ref] + alts
        qual = missing_to_none(lineparts[5])
        filt = missing_to_none(lineparts[6])
        info = {}
        for info_part in lineparts[7].split(';'):
            if '=' in info_part:
                key, value = info_part.split('=', 1)
                info[key] = value
            else:
                info[info_part] = None
        fmt = lineparts[8].split(':')

        data = [dict(zip(fmt, d.split(':'))) for d in lineparts[9:]]

        # Output format:
        # Sample   <line fields>    <alt allele fields>    <genotype fields>
        # Repeat as needed
        # Suppress lines if the relevant allele is not present in the given sample.

        for sample_i, sample in enumerate(samples):
            gts = [int(g) for g in data[sample_i]['GT'].split('/') if g != '.']
            for gt in gts:
                # Skip the reference GT.  This will naturally also lead to skipping
                # of 0/0 samples.
                if gt == 0:
                    continue

                if 'CSQ' in info:
                    vep_entries = info['CSQ'].split(',')
                    vep_parsed = [dict(zip(info_formats['CSQ'], v.split('|'))) for v in vep_entries]
                    vep = [v for v in vep_parsed if v['Allele'] == alts[gt - 1]]
                else:
                    vep = [{k: '' for k in info_formats['CSQ']}]

#                if 'Eigen' in info:
#                    eigen = [dict(zip(info_formats['Eigen'], info['Eigen'].split(',')[gt - 1].split('|')))]
#                else:
#                    eigen = [{}]
#                for eigen_field in info_formats['Eigen']:
#                    if eigen_field not in eigen[0]:
#                        eigen[0][eigen_field] = ''

#                if 'Cato3' in info:
#                    cato3 = [dict(zip(info_formats['Cato3'], info['Cato3'].split(',')[gt - 1].split('|')))]
#                else:
#                    cato3 = [{}]
#                for cato_field in info_formats['Cato3']:
#                    if cato_field not in cato3[0]:
#                        cato3[0][cato_field] = ''

#                if 'IARCTP53' in info:
#                    iarc = [dict(zip(info_formats['IARCTP53'], info['IARCTP53'].split(',')[gt - 1].split('|')))]
#                else:
#                    iarc = [{}]
#                for iarc_field in info_formats['IARCTP53']:
#                    if iarc_field not in iarc[0]:
#                        iarc[0][iarc_field] = ''

                for vep_entry in vep:
#                    for eigen_entry in eigen:
#                        for cato3_entry in cato3:
#                            for iarc_entry in iarc:
                                if 'AD' in data[sample_i]:
                                    ad = data[sample_i]['AD'].split(',')[gt]
                                else:
                                    ad = ''
#                                outstream.write('{s}\t{c}\t{p}\t{v}\t{r}\t{geno}\t{a}\t{ad}\t{dp}\t{ac}\t{af}'.format(
#                                    s = sample, c = chrom, p = pos, v = vid, r = ref, geno = '/'.join([alleles[i] for i in gts]), a = alts[gt - 1],
#                                    ad = ad, dp = data[sample_i]['DP'], 
#                                    ac = info['AC'].split(',')[gt - 1], af = info['AF'].split(',')[gt - 1]))
                                outstream.write('{s}\t{c}\t{p}\t{v}\t{r}\t{geno}\t{a}\t{ad}\t{dp}\t{tier}'.format(
                                    s = sample, c = chrom, p = pos, v = vid, r = ref, geno = '/'.join([alleles[i] for i in gts]), a = alts[gt - 1],
                                    ad = ad, dp = data[sample_i]['DP'], tier=info['tier']))

                                for vep_field in info_formats['CSQ']:
                                    outstream.write('\t{}'.format(vep_entry[vep_field]))
#                                for eigen_field in info_formats['Eigen']:
#                                    outstream.write('\t{}'.format(eigen_entry[eigen_field]))
#                                for cato_field in info_formats['Cato3']:
#                                    outstream.write('\t{}'.format(cato3_entry[cato_field]))
#                                for iarc_field in info_formats['IARCTP53']:
#                                    outstream.write('\t{}'.format(iarc_entry[iarc_field]))

                                outstream.write('\n')


if __name__ == '__main__':
    import sys

    samples, info_formats = parse_header(sys.stdin)

    process_data(sys.stdin, samples, info_formats, sys.stdout)
