import matplotlib.pyplot as plt
import numpy as np
import pod5 as p5
import sys
import struct
import os
import math
import pysam






# readToSignalPos = {}
# for line in open(sys.argv[1]):
#     line = line.rstrip().split('\t', 1)
#     readname = line[0]
#     mappings = line[1].split(',')
#     mappings = [[int(y) for y in x.split('.')] for x in mappings]
#     readToSignalPos[readname] = mappings

kmer_length = 9
readtoseq = {}
samfile = pysam.AlignmentFile(sys.argv[1], "rb")
c = 0
for s in samfile:
    if s.is_mapped and not s.is_supplementary and not s.is_secondary:
        alignstart = s.reference_start
        readname = s.query_name
        strand = -1 if s.is_reverse else 1
        chr = s.reference_name
        seq = s.query_sequence if not s.is_reverse else s.get_forward_sequence()
        # if s.is_reverse: seq = getrevcomp(seq)
        len_seq = len(seq) - kmer_length + 1  # to get the number of kmers

        ns = int(s.get_tag("ns")) ##number of signals
        ts = int(s.get_tag("ts")) ##signal start delay
        mv = s.get_tag("mv") ##signal list
        len_mv = len(mv)

        stride = mv[0]

        mvpos = 1
        move_count = 0

        start_signal_idx = ts
        currsiglen = 0
        # kmer_idx = 0

        kmersigpos = []

        ###Can we not save kmer_idx? Seems like each kmer will be represented in signal, so kmersigpos[i] = signal for kmer at pos i
        while mvpos < len_mv:
            mv_val = mv[mvpos]
            currsiglen += stride
            if mv_val == 1 or mv_val == len_mv-1:
                kmersigpos.append((start_signal_idx, start_signal_idx+currsiglen)) #(kmer_idx, start_signal_idx, start_signal_idx+currsiglen))
                start_signal_idx += currsiglen
                currsiglen = 0
                # kmer_idx += 1
            mvpos += 1
        readtoseq[readname] = [chr, strand, alignstart, seq, kmersigpos]

        c += 1
        if c % 1000 == 0: print('processed ' + str(c) + ' reads from .bam file')
print('done processing bam file')

signallines = []
readcodes = []
c = 0
out = open('.'.join(sys.argv[2].split('.')[:-1]) + '-scaled-kmersignal.bin', 'wb')

for pod5 in sys.argv[2:]:
    with p5.Reader(pod5) as reader:
        # read = next(reader.reads([selected_read_id]))
        for read in reader.reads():
            signal = read.signal_pa
            # if (not math.isnan(read.tracked_scaling.shift)) or (not math.isnan(read.tracked_scaling.scale)) or (not math.isnan(read.predicted_scaling.shift)) or (not math.isnan(read.predicted_scaling.scale)):
            #     print(read.tracked_scaling, read.predicted_scaling)
            # print(read.predicted_scaling.shift)

            scale = read.predicted_scaling.scale
            shift = read.predicted_scaling.shift
            signal = [(s-shift)/scale for s in signal]
            readname = str(read.read_id)
            if readname in readtoseq:
                kmerpos = 0
                for pos in readtoseq[readname][-1]:
                    signalLen = pos[1]-pos[0]
                    binary_data = struct.pack("i", c) + struct.pack("i", kmerpos) + struct.pack('i', signalLen) + struct.pack('f'*signalLen, *signal[pos[0]:pos[1]])
                    out.write(binary_data)
                    kmerpos += 1
                    # # possig = ','.join([str(x) for x in signal[pos[1]:pos[2]]])
                    # # out.write('\t'.join([str(c), str(pos[0]), possig]) + '\n')
                    # pickle.dump([c, pos[0], signal[pos[1]:pos[2]]], out)
                    # out.write('\t')
                    # signallines.append((c, pos[0], tuple(signal[pos[1]:pos[2]])))
                readcodes.append((readname, c))
                c += 1
                if c % 1000 == 0: print('processed ' + str(c) + ' reads from .pod5 file')
# out = open('.'.join(sys.argv[2].split('.')[:-1]) + '-kmersignal.pickle', 'wb')
# pickle.dump(signallines, out)
print('done processing pod5 files')
out.close()
out = open('.'.join(sys.argv[2].split('.')[:-1]) + '-scaled-readnameToCode.tsv', 'w')
for r in readcodes:
    out.write(str(r[1]) + '\t' + r[0] + '\n')
print('done writing readnames key')


