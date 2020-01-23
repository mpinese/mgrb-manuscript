#!/usr/bin/env python3

def ingest(path, name):
    print('\nATTACH "{path}" AS {name}_db;'.format(path=path, name=name))
    print('CREATE TABLE {name} AS SELECT chrom, pos, ref, alt, nRR as nRR_{name}, nRA as nRA_{name}, nAA as nAA_{name}, nmissing as nmissing_{name} FROM {name}_db.allelefreqs;'.format(name=name))
    print('DETACH {name}_db;'.format(name=name))
    print('CREATE UNIQUE INDEX {name}_idx ON {name}(chrom, pos, ref, alt);'.format(name=name))


def drop(name):
    print('\nDROP TABLE {name};'.format(name=name))


def join_and_clean(name_a, name_b, name_result):
    print('''
CREATE TABLE {name_result} AS
  SELECT A.chrom, A.pos, A.ref, A.alt, A.*, B.*
    FROM {name_a} A
    LEFT JOIN {name_b} B
    ON A.chrom = B.chrom AND A.pos = B.pos AND A.ref = B.ref AND A.alt = B.alt
  UNION ALL
  SELECT B.chrom, B.pos, B.ref, B.alt, A.*, B.*
    FROM {name_b} B
    LEFT JOIN {name_a} A
    ON A.chrom = B.chrom AND A.pos = B.pos AND A.ref = B.ref AND A.alt = B.alt
    WHERE A.chrom IS NULL;'''.format(name_a=name_a, name_b=name_b, name_result=name_result))
    print('DROP TABLE {name_a};'.format(name_a=name_a))
    print('DROP TABLE {name_b};'.format(name_b=name_b))
    print('CREATE UNIQUE INDEX {name_result}_idx ON {name_result}(chrom, pos, ref, alt);'.format(name_result=name_result))


def construct_join(inputs):
    # inputs is a list of (path, name) tuples
    assert len(inputs) > 1

    # Import the first input db.  We will progressively add to it.
    this_input_path, this_input_name = inputs[0]
    ingest(this_input_path, this_input_name)
    current_joined_table_name = this_input_name
    first_join = True

    # Progressively merge the input dbs with the growing joined db.
    for this_input_path, this_input_name in inputs[1:]:
        ingest(this_input_path, this_input_name)

        new_joined_table_name = current_joined_table_name + '_' + this_input_name
        join_and_clean(name_a=current_joined_table_name, name_b=this_input_name, name_result=new_joined_table_name)

        current_joined_table_name = new_joined_table_name

    # Construct the final select call.  This is what actually generates the output table.
    print('\nCREATE TABLE afwide AS SELECT chrom, pos, ref, alt,')
    for _, this_input_name in inputs[:-1]:
        print('  nRR_{this_input_name}, nRA_{this_input_name}, nAA_{this_input_name}, nmissing_{this_input_name},'.format(this_input_name=this_input_name))
    print('  nRR_{last_input_name}, nRA_{last_input_name}, nAA_{last_input_name}, nmissing_{last_input_name}'.format(last_input_name=inputs[-1][1]))
    print('FROM {new_joined_table_name};'.format(new_joined_table_name=new_joined_table_name))
    print('DROP TABLE {new_joined_table_name};'.format(new_joined_table_name=new_joined_table_name))
    print('\nCREATE UNIQUE INDEX afwide_idx on afwide(chrom, pos, ref, alt);')

#    print('\nCREATE TABLE aflong AS')
#    print('  SELECT chrom, pos, ref, alt, cohort, nRR, nRA, nAA, nmissing FROM afwide')
#    print('  CROSS APPLY ( VALUES')
#    for _, this_input_name in inputs[:-1]:
#        print('    ("{this_input_name}", nRR_{this_input_name}, nRA_{this_input_name}, nAA_{this_input_name}, nmissing_{this_input_name}),'.format(this_input_name=this_input_name))
#    print('    ("{last_input_name}", nRR_{last_input_name}, nRA_{last_input_name}, nAA_{last_input_name}, nmissing_{last_input_name})'.format(last_input_name=inputs[-1][1]))
#    print('  ) x (cohort, nRR, nRA, nAA, nmissing);')
#    print('\nCREATE UNIQUE INDEX aflong_idx ON aflong(chrom, pos, ref, alt, cohort);')

    print('\nVACUUM;')



if __name__ == '__main__':
    print('.echo on')
    construct_join([
       ('../MGRB.phase2final.GiaB_HCR.split.minrep.db', 'mgrborig'),
       ('../MGRB.phase2final.m.GiaB_HCR.split.minrep.db', 'mgrborig_m'),
       ('../MGRB.phase2final.f.GiaB_HCR.split.minrep.db', 'mgrborig_f'),
       ('../MGRB.phase2final.GiaB_HCR.split.minrep.somaticfiltered.db', 'mgrbsomfilt'),
       ('../MGRB.phase2final.ch10.GiaB_HCR.split.minrep.db', 'mgrborig_ch10'),
       ('../MGRB.phase2final.noch10.GiaB_HCR.split.minrep.db', 'mgrborig_noch10'),
       ('../MGRB.phase2final.ch10.GiaB_HCR.split.minrep.somaticfiltered.db', 'mgrbsomfilt_ch10'),
       ('../MGRB.phase2final.noch10.GiaB_HCR.split.minrep.somaticfiltered.db', 'mgrbsomfilt_noch10'),
       ('../UKBB.mf.less55.db', 'ukbb_0_55'),
       ('../UKBB.mf.55-60.db', 'ukbb_55_60'),
       ('../UKBB.mf.60-65.db', 'ukbb_60_65'),
       ('../UKBB.mf.65-70.db', 'ukbb_65_70'),
       ('../UKBB.mf.70-75.db', 'ukbb_70_75'),
       ('../UKBB.mf.more75.db', 'ukbb_75_inf'),
#       ('../ASRB.WGStier12.GiaB_HCR.split.minrep.db', 'asrborig'),
#       ('../ASRB.WGStier12.GiaB_HCR.split.minrep.m.db', 'asrborig_m'),
#       ('../ASRB.WGStier12.GiaB_HCR.split.minrep.f.db', 'asrborig_f'),
#       ('../ASRB.WGStier12.GiaB_HCR.split.minrep.somaticfiltered.db', 'asrbsomfilt'),
#       ('../ASRB.WGStier12.GiaB_HCR.split.minrep.somaticfiltered.m.db', 'asrbsomfilt_m'),
#       ('../ASRB.WGStier12.GiaB_HCR.split.minrep.somaticfiltered.f.db', 'asrbsomfilt_f'),
       ('../../../databases/GnomAD_hail/gnomad.genomes.r2.0.1.sites.autosomes.split.minrep.NFE.allelefreqs.db', 'gnomad'),
       ('../swegen_20180409.split.minrep.autosomes.allelefreqs.db', 'swegen'),
       ('../45andup.colorectalcancer.mf.GiaB_HCR.split.minrep.db', 'cancercrc'),
       ('../45andup.melanomacancer.mf.GiaB_HCR.split.minrep.db', 'cancermel'),
       ('../45andup.nonmelskincancer.mf.GiaB_HCR.split.minrep.db', 'cancernms'),
       ('../45andup.anycancer.mf.GiaB_HCR.split.minrep.db', 'cancerany'),
       ('../45andup.breastcancer.f.GiaB_HCR.split.minrep.db', 'cancerbrca_f'),
       ('../45andup.nocancer.mf.GiaB_HCR.split.minrep.db', 'nocancer'),
       ('../45andup.nocancer.f.GiaB_HCR.split.minrep.db', 'nocancer_f'),
       ('../45andup.nocancer.m.GiaB_HCR.split.minrep.db', 'nocancer_m'),
       ('../45andup.prostatecancer.m.GiaB_HCR.split.minrep.db', 'cancerpca_m')
    ])


