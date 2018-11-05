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

# 45 and up men with history of PCa
#target_samples = set([
#    'BAAAA','BAAAJ','BAAAM','BAAAP','BAAAX','BAABO','BAABT','BAABW',
#    'BAACE','BAACK','BAADM','BAAEB','BAAEH','BAAEV','BAAEZ','BAAFV',
#    'BAAFW','BAAGQ','BAAHD','BAAHM','BAAHR','BAAHY','BAAIK','BAAJJ',
#    'BAAJO','BAAKN','BAALK','BAAMJ','BAAND','BAANX','BAAOK','BAAOU',
#    'BAAOX','BAAOZ','BAAPV','BAAQC','BAAQW','BAAQY','BAAQZ','BAATW',
#    'BAAUA','BAAUC','BAAUI','BAAWF','BAAWW','BAAXB','BAAXD','BAAXQ',
#    'BAAYT','BAAZC','BAAZN','BAAZW','BABAE','BABAQ','BABAX','BABAZ',
#    'BABBB','BABBU','BABGL','BABGN','BABHF','BABHN','BABIV','BABJF',
#    'BABKH','BABKT','BABMV','BABNL','BABPD','BABPJ','BABPO','BABPQ',
#    'BABPX','BABQF','BABRE'])

# 45 and up men with no history of cancer
target_samples = set([
    'BAAAB','BAAAD','BAAAG','BAAAH','BAAAI','BAAAR','BAAAU','BAAAV',
    'BAABB','BAABC','BAABD','BAABL','BAABM','BAABN','BAABR','BAABS',
    'BAABU','BAABV','BAACG','BAACH','BAACI','BAACP','BAACU','BAACW',
    'BAACZ','BAADA','BAADB','BAADC','BAADE','BAADH','BAADI','BAADJ',
    'BAADK','BAADN','BAADP','BAADR','BAADS','BAADW','BAADX','BAADZ',
    'BAAEC','BAAEE','BAAEF','BAAEJ','BAAEK','BAAEP','BAAEU','BAAFA',
    'BAAFL','BAAFT','BAAFY','BAAGD','BAAGG','BAAGH','BAAGI','BAAGL',
    'BAAGO','BAAGV','BAAHB','BAAHE','BAAHF','BAAHH','BAAHI','BAAHJ',
    'BAAHL','BAAHO','BAAHQ','BAAHS','BAAHT','BAAHW','BAAHX','BAAIC',
    'BAAID','BAAIP','BAAIQ','BAAIS','BAAIU','BAAIX','BAAIY','BAAIZ',
    'BAAJA','BAAJD','BAAJE','BAAJF','BAAJL','BAAJM','BAAJQ','BAAJS',
    'BAAJV','BAAJW','BAAJY','BAAKE','BAAKJ','BAAKO','BAAKQ','BAAKU',
    'BAAKW','BAAKZ','BAALE','BAALH','BAALQ','BAALS','BAALU','BAALX',
    'BAALY','BAAMA','BAAMC','BAAMD','BAAMF','BAAML','BAAMN','BAAMQ',
    'BAAMS','BAAMT','BAAMU','BAAMV','BAAMX','BAAMZ','BAANF','BAANJ',
    'BAANK','BAANM','BAANO','BAANQ','BAANS','BAANT','BAANV','BAANZ',
    'BAAOA','BAAOC','BAAOG','BAAOJ','BAAOL','BAAON','BAAOO','BAAOW',
    'BAAOY','BAAPB','BAAPF','BAAPG','BAAPK','BAAPN','BAAPQ','BAAPR',
    'BAAPU','BAAPZ','BAAQD','BAAQO','BAAQS','BAAQT','BAAQV','BAARA',
    'BAARB','BAARC','BAARE','BAARF','BAARG','BAARI','BAARL','BAARS',
    'BAART','BAARY','BAASB','BAASC','BAASD','BAASJ','BAASM','BAASR',
    'BAASS','BAASW','BAASZ','BAATD','BAATF','BAATK','BAATP','BAATQ',
    'BAATS','BAATT','BAATU','BAATV','BAATX','BAAUD','BAAUE','BAAUJ',
    'BAAUK','BAAUL','BAAUQ','BAAUS','BAAUV','BAAVA','BAAVB','BAAVD',
    'BAAVG','BAAVI','BAAVK','BAAVL','BAAVM','BAAVS','BAAVT','BAAVW',
    'BAAVX','BAAVZ','BAAWA','BAAWH','BAAWP','BAAWQ','BAAWT','BAAWX',
    'BAAXI','BAAXJ','BAAXN','BAAXP','BAAXV','BAAXX','BAAXY','BAAYF',
    'BAAYJ','BAAYK','BAAYM','BAAYQ','BAAYR','BAAYS','BAAYU','BAAYW',
    'BAAZD','BAAZE','BAAZG','BAAZJ','BAAZL','BAAZM','BAAZO','BAAZR',
    'BAAZS','BAAZT','BAAZV','BABAA','BABAB','BABAH','BABAL','BABAM',
    'BABAO','BABAS','BABAT','BABAV','BABAW','BABBA','BABBE','BABBG',
    'BABBM','BABBN','BABBR','BABBS','BABBT','BABBV','BABBZ','BABCB',
    'BABCD','BABCE','BABCH','BABCI','BABCL','BABCM','BABCT','BABCV',
    'BABCW','BABCX','BABDM','BABDN','BABGF','BABGG','BABGR','BABGT',
    'BABGX','BABHB','BABHE','BABHK','BABHM','BABHO','BABHP','BABHQ',
    'BABHR','BABHV','BABHY','BABHZ','BABIC','BABID','BABIG','BABII',
    'BABIL','BABIM','BABIN','BABIR','BABIW','BABIZ','BABJE','BABJG',
    'BABJJ','BABJL','BABJR','BABJU','BABJV','BABJY','BABKI','BABKJ',
    'BABKO','BABKP','BABKQ','BABKS','BABKX','BABLA','BABLD','BABLE',
    'BABLF','BABLH','BABLI','BABLL','BABLO','BABLR','BABLT','BABLU',
    'BABLW','BABLX','BABMC','BABMD','BABME','BABMG','BABMI','BABMM',
    'BABMS','BABNA','BABNB','BABNC','BABND','BABNH','BABNS','BABNT',
    'BABNU','BABNV','BABNX','BABOA','BABOB','BABOH','BABOI','BABOK',
    'BABOL','BABOM','BABOQ','BABOR','BABOT','BABOV','BABOW','BABPE',
    'BABPF','BABPH','BABPM','BABPN','BABPS','BABPT','BABPU','BABPV',
    'BABPW','BABQA','BABQB','BABQD','BABQG','BABQN','BABQR','BABQX',
    'BABQY'])


# 45 and Up individuals with a history of CRC
#target_samples = set([
#    'BAAAE','BAAAM','BAAAV','BAABA','BAABJ','BAABN','BAACS','BAADP',
#    'BAADV','BAAGM','BAAGP','BAAGR','BAAGW','BAAHU','BAAHX','BAAIA',
#    'BAAIJ','BAAIQ','BAAKL','BAAKV','BAAMB','BAANW','BAAOJ','BAAQS',
#    'BAATH','BAAUR','BAAWO','BAAWX','BAAXD','BABAT','BABBE','BABBG',
#    'BABDQ','BABIW','BABLK','BABLQ','BABMG','BABMM','BABMR','BABNK',
#    'BABNR','BABNZ','BABOD','BABQD'])

# 45 and Up individuals with a history of melanoma
#target_samples = set([
#    'BAAAI','BAAAM','BAABK','BAACK','BAADN','BAAEA','BAAEZ','BAAFJ',
#    'BAAGB','BAAHB','BAAHG','BAAHS','BAAIC','BAAIM','BAAIP','BAAIQ',
#    'BAAJA','BAAJG','BAALQ','BAALV','BAAMS','BAANO','BAANS','BAAOW',
#    'BAAPI','BAAPW','BAAQY','BAARK','BAASD','BAAUE','BAAUV','BAAVF',
#    'BAAVJ','BAAVS','BAAXD','BAAXM','BAAXV','BAAYT','BAAZL','BABAM',
#    'BABAT','BABBM','BABBS','BABBX','BABCW','BABGR','BABJF','BABJR',
#    'BABLF','BABLI','BABNH','BABNX','BABOF','BABOG','BABOI','BABPM',
#    'BABPQ','BABQR','BABQW'])

# 45 and Up individuals with a history of non-melanoma skin cancer
#target_samples = set([
#    'BAACL','BAACS','BAADN','BAAGK','BAAGR','BAAGX','BAAHB','BAAHW',
#    'BAAIA','BAAIB','BAAIW','BAAJI','BAAJK','BAAKI','BAALT','BAAOZ',
#    'BAAPC','BAAPH','BAAPX','BAAQM','BAARD','BAARO','BAARV','BAATN',
#    'BAATP','BAAUI','BAAVT','BAAYZ','BABAK','BABBQ','BABBY','BABCD',
#    'BABDL','BABHM','BABMQ','BABNN','BABNX','BABOH','BABOO','BABOP',
#    'BABOW','BABPD','BABPL','BABPZ','BABQH','BABQP','BABQW','BABQX'])


gt_counter = collections.Counter()
i = 0
with open(inpath, 'rt') as infile:
    header = next(infile)
    sample_indices = [i for i, sid in enumerate(header.split('\t')[1:]) if sid in target_samples]
    for line in infile:
        vid, gts = line.rstrip().split('\t')
        chrom, pos, ref, alt = vid.split(':')
        gt_counter.clear()
        gt_counter.update([gts[i] for i in sample_indices])
        cur.execute('INSERT INTO allelefreqs VALUES (?, ?, ?, ?, ?, ?, ?, ?);', (chrom, pos, ref, alt, gt_counter['0'], gt_counter['1'], gt_counter['2'], gt_counter['.']))
        
        i += 1
        if (i == 1000000):
            i = 0
            print(vid)

conn.commit()
