#!/bin/bash

# This will convert BAM files from SRA
# from 10x data
# to fastq files that are compatible with cell ranger
# to actually be able to do the analysis

# first install bamtofastq
# https://github.com/10XGenomics/bamtofastq/releases
# and run it on your desired fastq
# what we did:

# get the fastq files from one GEM group (multiple libraries from this GEM group)
../bam/bamtofastq_macos ../bam/KD8_p1_0_possorted_genome_bam.bam.1 ../rawdata/