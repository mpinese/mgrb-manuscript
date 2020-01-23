#!/usr/bin/env python3
"""
Generate an allele frequency DB for MGRB variants used for PRS calculation.
"""
import sqlite3
import os.path
import logging


mgrb_dosages_path = '../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.popstar.dosages'
mgrb_afdb_path = '../MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.nocancer.over70.popstar.allele_freqs.dba'


logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO)

with open(mgrb_dosages_path, 'rt') as dosages_file:
    assert os.path.isfile(mgrb_afdb_path) == False, 'MGRB AF DB file {} already exists; refusing to overwrite.'.format(mgrb_afdb_path)
    conn = sqlite3.connect(mgrb_afdb_path)
    c = conn.cursor()
    c.execute('''CREATE TABLE mgrb (
        chrom     TEXT     NOT NULL,
        pos       INTEGER  NOT NULL,
        ref       TEXT     NOT NULL,
        alt       TEXT     NOT NULL,
        AC_MGRB   INTEGER  NOT NULL,
        AN_MGRB   INTEGER  NOT NULL,
        PRIMARY KEY (chrom, pos, ref, alt));''')

    sample_ids = next(dosages_file)

    last_chrom = None
    i = 0
    for line in dosages_file:
        vid, dosages = line.rstrip().split('\t', 1)
        chrom, pos, ref, alt = vid.split(':')

        if chrom != last_chrom:
            logging.info('Processing chromosome %s...', chrom)
            last_chrom = chrom
        i += 1
        if i % 100000 == 0:
            logging.info('%d variants processed', i)

        an = 0
        ac = 0
        for dosage in dosages:
            if dosage == '.':
                continue
            an += 2
            if dosage == '0':
                continue
            if dosage == '1':
                ac += 1
            elif dosage == '2':
                ac += 2
            else:
                assert False, 'Unexpected dosage character encountered'
        
        c.execute('INSERT INTO mgrb VALUES (?, ?, ?, ?, ?, ?);', (chrom, int(pos), ref, alt, ac, an))

    conn.commit()
    conn.close()

