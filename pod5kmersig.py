import matplotlib.pyplot as plt
import numpy as np
import pod5 as p5
import sys
import struct
import os


pod5 = sys.argv[1]

# readtoseq = {}
# last = None
# for line in open('can_reads.fasta'):
#     if line[0] == '>':
#         last = line[1:].rstrip()
#     else: readtoseq[last] = line.rstrip()


readToSignalPos = {}
for line in open(sys.argv[2]):
    line = line.rstrip().split('\t', 1)
    readname = line[0]
    mappings = line[1].split(',')
    mappings = [[int(y) for y in x.split('.')] for x in mappings]
    readToSignalPos[readname] = mappings

signallines = []
readcodes = []
c = 0
out = open('.'.join(sys.argv[2].split('.')[:-1]) + '-kmersignal.bin', 'wb')
with p5.Reader(pod5) as reader:
    # read = next(reader.reads([selected_read_id]))
    for read in reader.reads():
        signal = read.signal_pa
        readname = str(read.read_id)
        if readname in readToSignalPos:
            for pos in readToSignalPos[readname]:
                signalLen = pos[2]-pos[1]
                binary_data = struct.pack("i", c) + struct.pack("i", pos[0]) + struct.pack('i', signalLen) + struct.pack('f'*signalLen, *signal[pos[1]:pos[2]])
                out.write(binary_data)
                # # possig = ','.join([str(x) for x in signal[pos[1]:pos[2]]])
                # # out.write('\t'.join([str(c), str(pos[0]), possig]) + '\n')
                # pickle.dump([c, pos[0], signal[pos[1]:pos[2]]], out)
                # out.write('\t')
                # signallines.append((c, pos[0], tuple(signal[pos[1]:pos[2]])))
            readcodes.append((readname, c))
            c += 1
# out = open('.'.join(sys.argv[2].split('.')[:-1]) + '-kmersignal.pickle', 'wb')
# pickle.dump(signallines, out)
out.close()
out = open('.'.join(sys.argv[2].split('.')[:-1]) + '-readnameToCode.tsv', 'w')
for r in readcodes:
    out.write(str(r[1]) + '\t' + r[0] + '\n')



