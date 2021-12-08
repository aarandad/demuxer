# demuxer

This is how I run my pipeline

```
nextflow primerDemuxer_noerrors.nf
```

I modify the nextflow.config file with:
- a path to the folder containing fastq.gz files in reads
- the path to the fasta file with fwd or rev primers (either of the 3 subpools 1A, 1B or 2, one at a time) 
