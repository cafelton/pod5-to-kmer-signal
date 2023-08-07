import matplotlib.pyplot as plt
import sys, os, struct
import pysam
import numpy as np
from statistics import median, stdev



"""
python /private/groups/brookslab/cafelton/fusions-code/pod5kmersignal/plotSignalInRegion.py chr13:52310695-52310705
        can_mappings.reform-compressed-readnameToCode.tsv can_mappings.bam can_mappings.reform-compressed-kmersignal.bin 
        mod_mappings.reform-compressed-readnameToCode.tsv mod_mappings.bam mod_mappings.reform-compressed-kmersignal.bin 
"""

# python /private/groups/brookslab/cafelton/fusions-code/pod5kmersignal/plotSignalInRegion.py chr13:52310695-52310705 can_mappings.reform-compressed-readnameToCode.tsv can_mappings.bam can_mappings.reform-compressed-kmersignal.bin mod_mappings.reform-compressed-readnameToCode.tsv mod_mappings.bam mod_mappings.reform-compressed-kmersignal.bin

compbase = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
def getrevcomp(seq):
    seq = seq.upper()
    newseq = []
    for base in seq: newseq.append(compbase[base])
    return ''.join(newseq[::-1])


# readconvfile = sys.argv[2]
# bamfile = sys.argv[3]
# sigfile = sys.argv[4]

region = sys.argv[1]
# region = 'chr13:52310500-52310515' #'chr13:52310695-52310705' #'chr13:52308762-52308776'   #'chr13:52309895-52309915'
qstrand = 1
# modpos = 8

qchr, temp = region.split(':')
qstart, qend = temp.split('-')
qstart, qend = int(qstart), int(qend)

allfiles = []
for i in range(2, len(sys.argv), 3):
    allfiles.append([sys.argv[i], sys.argv[i+1], sys.argv[i+2]])


allreadtosig = []
for z in allfiles:#[[readconvfile, bamfile, sigfile], [sys.argv[5], sys.argv[6], sys.argv[7]]]:
    readconvfile, bamfile, sigfile = z[0], z[1], z[2]
    codetoread = {}
    for line in open(readconvfile):
        code, readname = line.rstrip().split('\t')
        code = int(code)
        codetoread[code] = readname

    readtoseq = {}
    regionrefseq = None
    samfile = pysam.AlignmentFile(bamfile, "rb")
    for s in samfile.fetch(qchr, qstart, qend):
        alignstart, alignend = s.reference_start, s.reference_end
        if alignstart <= qstart and alignend >= qend:
            readname = s.query_name
            strand = -1 if s.is_reverse else 1
            chr = s.reference_name
            seq = s.query_sequence#s.get_reference_sequence() #s.query_sequence ####UNSURE ABOUT THIS TBH
                ###Needs to be query sequence, but need to calculate query alignment to reference
            if not regionrefseq: regionrefseq = s.get_reference_sequence()[qstart-alignstart:qend-alignstart]
            if s.is_reverse: seq = getrevcomp(seq)

            ref, quer = alignstart, 0
            queryPosToRef = {}
            cigar = s.cigartuples
            for block in cigar:
                if block[0] in {0, 7, 8}:  # match, consumes both
                    for i in range(block[1]):
                        if qstart <= ref <= qend: queryPosToRef[quer] = ref
                        #if quer in ml: posOnGenome.append([ref + alignstart, ml[quer]])
                        ref += 1
                        quer += 1
                elif block[0] in {1, 4}:  # consumes query
                    quer += block[1]
                elif block[0] in {2, 3}:  # consumes reference
                    ref += block[1]

            readtoseq[readname] = [chr, strand, alignstart, seq, queryPosToRef]


    filesize = os.stat(sigfile).st_size
    readtosiginregion = {}
    c = 0
    with open(sigfile, mode='rb') as file:
        while c < filesize:
            rcode, kmeridx, signalcount = struct.unpack('iii', file.read(12))
            readname = codetoread[rcode]
            if readname in readtoseq:
                chr, strand, alignstart, seq, queryPosToRef = readtoseq[readname]
                # print(qchr, qstart, qend, qstrand, '\t', chr, strand, kmeridx+alignstart)
                if qchr == chr and qstrand == strand and kmeridx in queryPosToRef: #qstart <= kmeridx + alignstart <= qend:
                    if readname not in readtosiginregion: readtosiginregion[readname] = {}#{'sig':[], 'seq': ''}
                    sig = struct.unpack('f' * signalcount, file.read(signalcount * 4))
                    readtosiginregion[readname][queryPosToRef[kmeridx]] = [seq[kmeridx]] + list(sig)
                    # readtosiginregion[readname]['sig'].append(sig)
                    # readtosiginregion[readname]['seq'] += seq[kmeridx]

                else: file.read(signalcount*4)
            else: file.read(signalcount * 4)
            c += 12 + signalcount * 4
    allreadtosig.append(readtosiginregion)

print('files processed, plotting')

###non-overlapping plots
# fig, axs = plt.subplots(len(readtosiginregion))
# c = 0
# # print(readtosiginregion)
# for r in readtosiginregion:
#     linepos = 0
#     for p in range(qstart, qend):
#         if p in readtosiginregion[r]:
#             # print(r, p, readtosiginregion[r][p][0], len(readtosiginregion[r][p])-1)
#             siglen = len(readtosiginregion[r][p])-1
#             # axs[c].axvline(linepos + siglen, c='black', linewidth=1)
#             axs[c].text(linepos + siglen / 2, 100, readtosiginregion[r][p][0], ha='center')
#             axs[c].plot(list(range(linepos, siglen+linepos)), readtosiginregion[r][p][1:])
#         else:
#             axs[c].plot(list(range(linepos, 40+linepos)), [70]*40)
#         linepos += 40#siglen
#     c += 1


fig = plt.figure()
d = 0
sigdiff = [[] for x in range(qend-qstart)]
for readtosiginregion in allreadtosig:
    c = 0
    avgsig = [[] for x in range(qend-qstart)]
    for r in readtosiginregion:
        linepos = 0
        for p in range(qstart, qend):
            if p in readtosiginregion[r] and regionrefseq[p-qstart] == readtosiginregion[r][p][0].upper():
                # print(r, p, readtosiginregion[r][p][0], len(readtosiginregion[r][p])-1)
                siglen = len(readtosiginregion[r][p])-1
                # axs[c].axvline(linepos + siglen, c='black', linewidth=1)
                xaxis = list(np.linspace(linepos+qstart, linepos + 1 + qstart, siglen + 1))[:-1]
                # axs[c].text(linepos + siglen / 2, 100, readtosiginregion[r][p][0], ha='center')
                plt.plot(xaxis, readtosiginregion[r][p][1:], c = 'C' + str(d), alpha=0.2)
                if c == 0:
                    plt.text(qstart + linepos + .5, 100, regionrefseq[p-qstart], ha='center')
                avgsig[p-qstart].append(sum(readtosiginregion[r][p][1:])/siglen)
            linepos += 1#siglen
        c += 1
    # for p in range(qend-qstart):
    #     plt.plot([p+qstart, qstart+p+1], [median(avgsig[p])]*2, c='C' + str(d), alpha=1)
    #     # plt.text(qstart + p + .5, 50 + d * 80 + ((p%2)*5), str(round(stdev(avgsig[p]))), ha='center', c='C' + str(d))
    #     sigdiff[p].append(stdev(avgsig[p]))
    d += 1
# if len(sys.argv) >= 6:
#     for p in range(qend-qstart):
#         plt.text(qstart + p + .5, 130, str(round(abs(sigdiff[p][1]-sigdiff[p][0]))), ha='center')

plt.savefig(bamfile.split('.')[0] + region + '-readsigplots-overlaid-with-variability.png', dpi=600)
print('done')
