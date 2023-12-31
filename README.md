# pod5-to-kmer-signal
Tools for getting signals aligned to read sequence for R10 flow cells and Dorado



You need to start with a merged .pod5 current signal file and a basecalled and sorted and indexed .bam file that contains both the MD tag and the move table (Dorado --emit_moves).

This tool was inspired by squigalizer: https://github.com/hiruna72/squigualiser

# NEW VERSION (faster, more compressed)
```bash
usage: python[3+] bampod5kmersig-arrow.py -b sample.bam -p sample.pod5 -n 2000 -c chrI
Predicting the best modification threshold and dist from positive and negative control bam files. Will output a tsv file containing those values.
options:
  -h, --help            show this help message and exit
  -b B, --bam B         R10 basecalled bam file, should have moves table and be indexed. Can be filtered to locus or not. Only reads in this bam will have signal called for kmers
  -p P, --pod5 P        R10 pod5 file, should have all reads that are included in bam. This should be a single file, not multiple
  -o O, --outputprefix O
                        Output filename prefix
  -s, --scale           whether to scale the signal values to facillitate comparison between reads, requires that pod5 reads have a non-nan value in the predicted_scaling area
  -c C, --chr C         chromosome to get reads from
  -n N, --numreads N    number of reads to get. if -c is specified, will be first n reads on the chromosome specified, otherwise will be first n reads on first chr
```

In the output there will be one row for each base in each read.

The schema for the parquet table is as follows:
```bash
readcode: string
kmeridx: int64
qkmer: string (9 chars long)
refpos: int64
refkmer: string (1 char)
signalLen: int64
signal: list<item: double>
  child 0, item: double
```



To look at a specific region in the output:

look at current signal for one file at a specific locus:
```bash
python3 plotSignalInRegion-arrow.py chrI:4300-4320 sample1.bam sample1-kmersignal-complete.parquet
```

look at overlapping current signal for multiple files at a specific locus:
```bash
python3 plotSignalInRegion-arrow.py chrI:4300-4320 sample1.bam sample1-kmersignal-complete.parquet sample2.bam sample2-kmersignal-complete.parquet
```





# OLD INSTRUCTIONS FOR OLD VERSION

I reccommend installing squigalizer with the conda environment:
```bash
git clone https://github.com/hiruna72/squigualiser.git
cd squigualiser
conda create -n squig python=3.8.0 -y
conda activate squig
python setup.py install
squigualiser --help
```

You can then start by extracting the move data from the bam file by running:
```bash
squigualiser reform --sig_move_offset 0 --kmer_length 9 --bam mappings.bam -o mappings.reform.tsv
```

This will output a fairly large file, which can then be compressed with:
```bash
python compressSquigReform.py mappings.reform.tsv
rm mappings.reform.tsv
```

You can then extract the signal data from the pod5 files with:
```bash
python pod5kmersig.py mappings.reform-compressed.tsv reads.pod5
rm mappings.reform-compressed.tsv
```

If your pod5s are in multiple files, you can run:
```bash
python pod5kmersig.py mappings.reform-compressed.tsv *.pod5
```

You will now have two files: -readnameToCode.tsv and -kmersignal.bin. To extract all of the data that is available through eventalign from these files, you can use the file processpod5kmersig.py as your template. Each kmer in each read can be processed and analyzed individually if you want.

If you want to visualize your signal in a genomic region, you can run the below code. This does force the dwell time in the pore for each base in each read to be the same on the plot, so if you want to look at more raw data, go into the file and uncomment the section labeled ###non onverlapping plots while commenting out the ###overlaid plots section.
```bash
python /private/groups/brookslab/cafelton/fusions-code/pod5kmersignal/plotSignalInRegion.py chr13:52309895-52309915 can_mappings.reform-compressed-readnameToCode.tsv can_mappings.bam can_mappings.reform-compressed-kmersignal.bin 
```

You can also plot multiple files on the same axis - they will plot in different colors.
```bash
python /private/groups/brookslab/cafelton/fusions-code/pod5kmersignal/plotSignalInRegion.py chr13:52309895-52309915 can_mappings.reform-compressed-readnameToCode.tsv can_mappings.bam can_mappings.reform-compressed-kmersignal.bin mod_mappings.reform-compressed-readnameToCode.tsv mod_mappings.bam mod_mappings.reform-compressed-kmersignal.bin
```
