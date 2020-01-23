cancer_samples = set([x.strip() for x in open('../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/45andup_followup_cancer.sample_list', 'rt')])
young_samples = set([x.strip() for x in open('../../../GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/45andup_followup_under70.sample_list', 'rt')])

infile = open('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.popstar.dosages', 'rt')
outfile_mgrb = open('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.popstar.dosages', 'wt')
outfile_cancer = open('../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.cancer.over70.popstar.dosages', 'wt')

header = next(infile).rstrip().split('\t')
dosage_samples = set(header[1:])
cancer_samples = cancer_samples & dosage_samples
young_samples = young_samples & dosage_samples
mgrb_samples = dosage_samples - cancer_samples - young_samples
print(len(mgrb_samples))
cancer_samples = list(cancer_samples)
mgrb_samples = list(mgrb_samples)
young_samples = list(young_samples)
cancer_samples.sort()
mgrb_samples.sort()
young_samples.sort()
young_indices = [header.index(sampleid) for sampleid in young_samples]
cancer_indices = [header.index(sampleid) for sampleid in cancer_samples]
mgrb_indices = [header.index(sampleid) for sampleid in mgrb_samples]

outfile_mgrb.write('vid\t' + '\t'.join([header[idx] for idx in mgrb_indices]) + '\n')
outfile_cancer.write('vid\t' + '\t'.join([header[idx] for idx in cancer_indices]) + '\n')

for line in infile:
    variant, dosages = line.rstrip().split('\t')
    outfile_mgrb.write(variant + '\t' + ''.join([dosages[idx-1] for idx in mgrb_indices]) + '\n')
    outfile_cancer.write(variant + '\t' + ''.join([dosages[idx-1] for idx in cancer_indices]) + '\n')

