# Check-In 1
## Addressing Prior Feedback
- datasets: 
    - To recreate a portion of the raw data processing we will use the following subset offastq files:
        - fastq file name,SRA Run ID
        - KD8_seq1_p1_mRNA_0_S1_L004,SRR19330988
        - KD8_seq1_p1_sgRNA_0_S1_L004,SRR19330979
        - KD8_seq2_p1_mRNA_0_S1_L004,SRR19330980
        - KD8_seq2_p1_sgRNA_0_S1_L004,SRR19330971 
    - To recreate the processing of aligned files to cell-gene read count matrices, we will use the follow BAM file:
        - KD8_p1_0_possorted_genome_bam.bam.1
- TAPseq: 
    - TAPseq will allow us access off target effects to identify strongly knocked down genes based on the identity of the gRNA. The functionality we plan to implement is to write code that looks at the knocked down genes, and flag any genes that are within 2 kb of the gRNA targeting site, and use these as candidate off-target genes.
- Goals: 
    - we plan to recreate figures 1B, 2B, 2C, and maybe 3B
## Progress
Dan Le:
- investigated raw fastq files and attempted to recreate 10x cell-ranger pipeline on subset of raw sequencing data
- pivoted to processing of bam/sam files when it turned out that SRA fastq files aren't compatible with cell ranger
- dumped sam files to laptop, filtered UMIs that pass QC check, and get read counts for each gene per cell
Jander Kugelman:
- rewrote a good portion of perturbseq from python to R
- rewritten (not validated):
    - cell_cycle.R
    - differential_expression.R
    - expression_normalization.R
    - util.R
    - zzz.R
- skeleton:
    - cell_population.R
    - transformers.R
## Project Organization
- check in writups in /doc
- perturbseq rewrite in /src/perturbseqR
    - rewritten R scripts in /src/perturbseqR/R
## Struggles/Questions
- cell-ranger does not accept sra fastq files (no flow cell ID)
    - reason for switching to bam/sam processing