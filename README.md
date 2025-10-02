# Transcriptome Correlation and Off-Target Analysis in Perturb-seq

## Description (50-100 words)


## Example Published Figure: S2, Replogle et al. 2022

![figS2](src/sup_figs/figS2.jpg)

## Dataset IDs:
FASTQ and BAM:

Raw Data: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA831566

Metadata: https://plus.figshare.com/articles/dataset/_Mapping_information-rich_genotype-phenotype_landscapes_with_genome-scale_Perturb-seq_Replogle_et_al_2022_SRA_and_GEO_file_manifest/20022944 

Database for stretch goal: https://thebiogrid.org/

## Software:
- Python 2.7
    - PerturbSeq Library: https://github.com/thomasmaxwellnorman/Perturbseq_GI/tree/master
- R 4.5.0
    - write PerturbSeq library in R
- CellRanger 4.0.0 : https://www.10xgenomics.com/support

## Goals:

1. Recreate data-containing figures/panels from Replogle et al. 2022 (15 figures: 7 main text, 8 supplemental) using R
    - Code from paper is developed in python
2. Transcriptome mapping: correlate distance between genes and determine off target KD potential of paired genes
    - Start with bidirectional reporters and look for off-target KD
    - See if KD genes have correlated transcriptome readouts
    - Assign function to unannotated genes
3. (Stretch) Use raw data from CRISPR/CRISPRi screens and pipeline generated for figure recreation to map genotype-phenotype landscapes in astrocytes with PerturbSeq
    - Potentially follow-up with transcriptome mapping

## References:
Replogle, Joseph M., Reuben A. Saunders, Angela N. Pogson, et al. 2022. “Mapping Information-Rich Genotype-Phenotype Landscapes with Genome-Scale Perturb-Seq.” Cell 185 (14): 2559-2575.e28. https://doi.org/10.1016/j.cell.2022.05.013.
