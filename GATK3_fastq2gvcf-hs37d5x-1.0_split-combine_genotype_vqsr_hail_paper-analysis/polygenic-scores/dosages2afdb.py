import collections
import sqlite3
import sys

inpath, outpath = sys.argv[1:]

conn = sqlite3.connect(outpath)
cur = conn.cursor()
cur.execute('''CREATE TABLE allelefreqs (
    chrom    TEXT    NOT NULL,
    pos      INTEGER NOT NULL,
    ref      TEXT    NOT NULL,
    alt      TEXT    NOT NULL,
    nRR      INTEGER NOT NULL,
    nRA      INTEGER NOT NULL,
    nAA      INTEGER NOT NULL,
    nmissing INTEGER NOT NULL,
    PRIMARY KEY (chrom, pos, ref, alt));''')

gt_counter = collections.Counter()
i = 0
with open(inpath, 'rt') as infile:
    header = next(infile)
    for line in infile:
        vid, gts = line.rstrip().split('\t')
        chrom, pos, ref, alt = vid.split(':')
        gt_counter.clear()
        gt_counter.update(gts)
        cur.execute('INSERT INTO allelefreqs VALUES (?, ?, ?, ?, ?, ?, ?, ?);', (chrom, pos, ref, alt, gt_counter['0'], gt_counter['1'], gt_counter['2'], gt_counter['.']))
        
        i += 1
        if (i == 1000000):
            i = 0
            print(vid)

conn.commit()
