import matplotlib.pyplot as plt
import numpy as np
import pod5 as p5

# Using the example pod5 file provided
pod5 = "can_reads.pod5"
selected_read_id = '6e37823a-9398-4be8-b111-65cab029f4e0'

readtoseq = {}
last = None
for line in open('can_reads.fasta'):
    if line[0] == '>':
        last = line[1:].rstrip()
    else: readtoseq[last] = line.rstrip()

print(len(readtoseq[selected_read_id]))

plt.figure()
for line in open('can_mappings_kmerlen1.reform-compressed.tsv'):
    line = line.rstrip().split('\t', 1)
    if line[0] == selected_read_id:
        signals = line[1].split(',')[:100]
        signals = [[int(y) for y in x.split('.')] for x in signals]
        # print(signals)
        startline = 0
        c = 0
        lasttexthigh = False
        while startline <= 400:
            if signals[c][1] >= 200:
                plt.axvline(signals[c][1], c='black', linewidth=1)
                plt.axvline(signals[c][2], c='black', linewidth=1)
                textheight = 115 if lasttexthigh else 125
                lasttexthigh = False if lasttexthigh else True
                plt.text((signals[c][1] + signals[c][2])/2, textheight, readtoseq[selected_read_id][signals[c][0]], ha='center')
            c += 1
            startline = signals[c][1]


with p5.Reader(pod5) as reader:

    # Read the selected read from the pod5 file
    # next() is required here as Reader.reads() returns a Generator
    read = next(reader.reads([selected_read_id]))

    # Get the signal data and sample rate
    sample_rate = read.run_info.sample_rate
    signal = read.signal_pa

    # Compute the time steps over the sampling period
    time = range(len(signal))#np.arange(len(signal)) / sample_rate
    print(len(signal))
    # Plot using matplotlib
    print(signal[:10])
    plt.plot(time[200:400], signal[200:400])
plt.savefig(pod5.split('.')[0] + '_' + selected_read_id + '_200-400sig_withbasesep.png', dpi=600)