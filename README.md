# pod5-to-kmer-signal
Tools for getting signals aligned to read sequence for R10 flow cells and Dorado



You need to start with .pod5 current signal files and a basecalled .bam file that contains both the MD tag and the move table (Dorado --emit_moves).

This tool depends on and was inspired by squigalizer: https://github.com/hiruna72/squigualiser

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
