import gzip
import warnings
import pprint
import pysam


def variant_clinvar_to_vcf(ref_clinvar, alt_clinvar, chrom, start_clinvar, stop_clinvar, reference):
    """Convert a ClinVar-format variant to VCF.  ref_clinvar, alt_clinvar, and chrom
    are strings, start_clinvar and stop_clinvar are integers, and reference is a
    pysam.FastaFile object.  Returns a tuple of (pos, ref, alt) on success, or
    (None, None, None) on failure.

    Note that we don't go to extreme lengths to ensure standards-compliance, in
    particular there is no attempt to justify or minimise alleles.  It is assumed
    that a subsequent vt norm stage will be run.
    """
    
    # The ClinVar table format lists reference and alternate alleles
    # (fields ReferenceAllele, AlternateAllele).  In most cases these
    # can be directly translated to the VCF alleles, but in the case of
    # simple indels the ClinVar encoding is '-', which is not VCF-
    # compatible.  In those cases we need to source the 'missing'
    # base from the reference file.

    if is_simple_var(ref_clinvar, alt_clinvar):
        # Simple variant
        return simple_var_clinvar_to_vcf(ref_clinvar, alt_clinvar, chrom, start_clinvar, stop_clinvar, reference)
    elif is_simple_ins(ref_clinvar, alt_clinvar):
        # Simple insertion
        return simple_ins_clinvar_to_vcf(ref_clinvar, alt_clinvar, chrom, start_clinvar, stop_clinvar, reference)
    elif is_simple_del(ref_clinvar, alt_clinvar):
        # Simple deletion
        return simple_del_clinvar_to_vcf(ref_clinvar, alt_clinvar, chrom, start_clinvar, stop_clinvar, reference)
    else:
        # Conversion failed -- unknown variant type
        return (None, None, None)


def is_valid_vcf_allele(allele_string):
    for base in allele_string:
        if (base != 'G' and base != 'C' and base != 'A' and base != 'T'):
            return False
    return True


def is_simple_var(ref_clinvar, alt_clinvar):
    return is_valid_vcf_allele(ref_clinvar) and is_valid_vcf_allele(alt_clinvar)


def is_simple_ins(ref_clinvar, alt_clinvar):
    return ref_clinvar == '-' and is_valid_vcf_allele(alt_clinvar)


def is_simple_del(ref_clinvar, alt_clinvar):
    return is_valid_vcf_allele(ref_clinvar) and alt_clinvar == '-'


def simple_var_clinvar_to_vcf(ref_clinvar, alt_clinvar, chrom, start_clinvar, stop_clinvar, reference):
    fasta_ref = reference.fetch(reference=chrom, start=start_clinvar-1, end=stop_clinvar)
    if ref_clinvar != fasta_ref:
        warnings.warn('Assertion failed: reference sequence mismatch, for variant R:{} A:{} C:{} P:{}-{} F:{}'.format(ref_clinvar, alt_clinvar, chrom, start_clinvar, stop_clinvar, fasta_ref))
        return (None, None, None)
    return (start_clinvar, ref_clinvar, alt_clinvar)


def simple_ins_clinvar_to_vcf(ref_clinvar, alt_clinvar, chrom, start_clinvar, stop_clinvar, reference):
    # The insert is *between* start_clinvar and stop_clinvar, in 1-based coordinates.
    # This is similar to HGVS representation.
    # To convert to VCF, fetch the reference base immediately upstream of the insert
    # site (ie at start_clinvar), and prepend it to both alleles.  The VCF variant location
    # is then start_clinvar.
    if stop_clinvar != start_clinvar + 1:
        warnings.warn('Assertion failed: stop != start + 1, for variant R:{} A:{} C:{} P:{}-{}'.format(ref_clinvar, alt_clinvar, chrom, start_clinvar, stop_clinvar))
        return (None, None, None)
    upstream_ref = reference.fetch(reference=chrom, start=start_clinvar-1, end=stop_clinvar-1)
    return (start_clinvar, upstream_ref, upstream_ref + alt_clinvar)


def simple_del_clinvar_to_vcf(ref_clinvar, alt_clinvar, chrom, start_clinvar, stop_clinvar, reference):
    # The deletion spans start_clinvar to stop_clinvar, in 1-based inclusive coordinates.
    # To convert to VCF, fetch the reference base immediately upstream of the first
    # deleted base (ie at start_clinvar - 1, in 1-based coordinates), and prepend it to
    # both alleles.  The VCF variant location is then start_clinvar - 1.
    if stop_clinvar != start_clinvar + len(ref_clinvar) - 1:
        warnings.warn('Assertion failed: stop != start + len(ref) - 1, for variant R:{} A:{} C:{} P:{}-{}'.format(ref_clinvar, alt_clinvar, chrom, start_clinvar, stop_clinvar))
        return (None, None, None)
    fasta_ref = reference.fetch(reference=chrom, start=start_clinvar-1, end=stop_clinvar)
    if ref_clinvar != fasta_ref:
        warnings.warn('Assertion failed: reference sequence mismatch, for variant R:{} A:{} C:{} P:{}-{} F:{}'.format(ref_clinvar, alt_clinvar, chrom, start_clinvar, stop_clinvar, fasta_ref))
        return (None, None, None)
    upstream_ref = reference.fetch(reference=chrom, start=start_clinvar-2, end=start_clinvar-1)
    return (start_clinvar - 1, upstream_ref + ref_clinvar, upstream_ref)


def convert_clinvar_clinsig(clinsig_string):
    return ','.join((clinsig_part.title().replace(' ', '') for clinsig_part in clinsig_string.split(',')))


def convert_clinvar_reviewstatus(reviewstatus_string):
    return reviewstatus_string.replace(',', '').title().replace(' ', '')


def convert_vcf(infile, outfile, ref):
    header = next(infile).rstrip().split('\t')
    assert header[0].startswith('#')
    header[0] = header[0][1:]

    outfile.write('##fileformat=VCFv4.2\n')
    outfile.write('##INFO=<ID=AlleleID,Number=1,Type=Integer,Description="integer value as stored in the AlleleID field in ClinVar">\n')
    outfile.write('##INFO=<ID=ClinicalSignificance,Number=.,Type=String,Description="clinical significance reported for this variant">\n')
    outfile.write('##INFO=<ID=ClinSigSimple,Number=1,Type=Integer,Description="0 = no current value of Likely pathogenic/Pathogenic, 1 = at least one current submissions interpreting as Likely pathogenic/Pathogenic">\n')
    outfile.write('##INFO=<ID=rsid,Number=1,Type=Integer,Description="rsID in dbSNP">\n')
    outfile.write('##INFO=<ID=ReviewStatus,Number=1,Type=String,Description="highest review status for reporting this measure">\n')
    outfile.write('##INFO=<ID=NumberSubmitters,Number=1,Type=Integer,Description="number of submitters describing this variant">\n')
    outfile.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

    converted = 0
    skipped = 0
    failed = 0

    for line in infile:
        fields = dict(zip(header, line.rstrip().split('\t')))

        if fields['Assembly'] != 'GRCh37':
            continue

        if fields['ReferenceAllele'] == fields['AlternateAllele']:
            skipped += 1
            continue

        if fields['ReferenceAllele'] == 'na' or fields['AlternateAllele'] == 'na':
            skipped += 1
            continue

        if fields['ClinicalSignificance'] == '-':
            skipped += 1
            continue

        clinvar_start = int(fields['Start'])
        clinvar_stop = int(fields['Stop'])
  
        vcf_variant = variant_clinvar_to_vcf(fields['ReferenceAllele'], fields['AlternateAllele'], fields['Chromosome'], clinvar_start, clinvar_stop, ref)

        if vcf_variant[0] == None:
            warnings.warn('Could not convert variant: rs{} (ref: {}, alt: {})'.format(fields['RS# (dbSNP)'], fields['ReferenceAllele'], fields['AlternateAllele']))
            failed += 1
            continue

        outfile.write('{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filt}\t'.format(
            chrom=fields['Chromosome'], pos=vcf_variant[0], id='.', 
            ref=vcf_variant[1], alt=vcf_variant[2], qual='.', filt='.'))
        outfile.write('AlleleID={aid};ClinicalSignificance={cs};ClinSigSimple={css};rsid={rsid};ReviewStatus={revstat};NumberSubmitters={numsub}\n'.format(
            aid=fields['AlleleID'], cs=convert_clinvar_clinsig(fields['ClinicalSignificance']), css=fields['ClinSigSimple'], 
            rsid=fields['RS# (dbSNP)'], revstat=convert_clinvar_reviewstatus(fields['ReviewStatus']), numsub=fields['NumberSubmitters']))
        converted += 1
    warnings.warn('Converted: {}   Skipped: {}   Failed: {}'.format(converted, skipped, failed))


if __name__ == '__main__':
    import sys

    # ../source_data/variant_summary.txt.gz
    infile = gzip.open(sys.argv[1])
    # ../../../resources/hs37d5x/reference_genome/hs37d5x.fa
    ref = pysam.FastaFile(sys.argv[2])

    convert_vcf(infile, sys.stdout, ref)

