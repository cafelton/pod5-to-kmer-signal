import sys, os, struct
import pysam

readconvfile = sys.argv[1]
fastafile = sys.argv[2]
sigfile = sys.argv[3]

codetoread = {}
for line in open(readconvfile):
    code, readname = line.rstrip().split('\t')
    codetoread[code] = readname

# readtoseq = {}
# last = None
# for line in open(fastafile):
#     if line[0] == '>':
#         last = line[1:].rstrip()
#     else: readtoseq[last] = line.rstrip()

compbase = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
def getrevcomp(seq):
    newseq = []
    for base in seq: newseq.append(compbase[base])
    return ''.join(newseq[::-1])

readtoseq = {}
samfile = pysam.AlignmentFile(args.b, "rb")
for s in samfile:
    alignstart = s.reference_start
    readname = s.query_name
    strand = -1 if s.is_reverse else 1
    chr = s.reference_name
    seq = s.query_sequence
    if s.is_reverse: seq = getrevcomp(seq)
    readtoseq[readname] = [chr, strand, alignstart, seq]

readcode, siglist, siglenlist = None, [], []
filesize = os.stat(sigfile).st_size
c = 0
with open(sigfile, mode='rb') as file:
    while c < filesize:
        rcode, kmeridx, signalcount = struct.unpack('iii', file.read(12))
        sig = struct.unpack('f'*signalcount, file.read(signalcount*4))
        c += 12 + signalcount * 4
        if rcode != readcode:
            if readcode:
                readname = codetoread[readcode]
                chr, strand, alignstart, seq = readtoseq[readname]
                ### This is where you could pass the data to the rest of your functions
                """
                metadata, scores = modPredict(read=readname, chrom=chr, start=alignstart,
                                              seq = seq, signals = siglist,
                                              siganlLengthList=siglenlist, strandInt=strand, 
                                              downScores = downScores, method = method)
                """
            readcode = rcode
            siglist, siglenlist = [], []
        siglist.extend(sig)
        siglenlist.append(len(sig))


# for line in open(sigfile):
#     line = line.rstrip().split('\t')
#     sig = [float(x) for x in line[2].split(',')]
#     rcode = int(line[0])
#     if rcode != readcode:
#         if readcode:
#             readname = codetoread[readcode]
#             chr, strand, alignstart, seq = readtoseq[readname]
#             ### This is where you could pass the data to the rest of your functions
#             """
#             metadata, scores = modPredict(read=readname, chrom=chr, start=alignstart,
#                                           seq = seq, signals = siglist,
#                                           siganlLengthList=siglenlist, strandInt=strand,
#                                           downScores = downScores, method = method)
#             """
#         readcode = line[0]
#         siglist, siglenlist = [], []
#     siglist.extend(sig)
#     siglenlist.append(len(sig))