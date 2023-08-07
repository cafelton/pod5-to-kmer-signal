import sys


readname = None
pos = []
out = open('.'.join(sys.argv[1].split('.')[:-1]) + '-compressed.tsv', 'w')
for line in open(sys.argv[1]):
    line = line.rstrip().split('\t')
    if line[0] != readname:
        if readname:
            out.write(readname + '\t' + ','.join(pos) + '\n')
        readname = line[0]
        pos = []
    pos.append('.'.join(line[1:]))
out.close()