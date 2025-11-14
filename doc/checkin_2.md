---
editor_options: 
  markdown: 
    wrap: 72
---

# Check-In 2

## Addressing Prior Feedback

-   How will you test and document tests in this repo about the
    rewritting of PerturbSeq?

-   The original PerturbSeq documentation contains ipynb notebooks with
    sample usage of the package. Additionally, there is a separate repo
    that is self-contained and has all the data required for a
    perturbseq demo, which we think is the best initial test:
    <https://github.com/thomasmaxwellnorman/perturbseq_demo> (i.e., we
    will write the perturbseq_demo Python notebook into an R notebook
    and test to see if the demo can be recreated)

## Progress 

Dan Le:

- To recreate the processing of raw sequecing data on a subset of data, we obtained raw FASTQ files using bamtofastq (<https://github.com/10XGenomics/bamtofastq>). With this, we managed to convert the bam file associated with library KD8_p1_0 (K562 8 Day Genome Wide Peturb-seq screen) to the original fastq files, generating all the raw sequencing reads for gene transcripts as well as the captured guide RNAs for approximately 15,000 cells. Uploading this to 10x's cloud CellRanger, we managed to align to the hg38 as done in Replogle et al. 2022, retrieving gene expression data from almost 15000 unique cells. A summary of a pilot alignment is stored in /processed/cell_ranger_output/20251107/. Unfortunately, the authors did not provide a Cell Ranger compatible Feature Reference CSV for CRISPR guide capture, which provides the CRISPR sgRNA sequences used to assign the perturbed gene to each cell based on the reads from the CRISPR guide capture library. I have obtained the protospacer sequences for the dual-sgRNA vectors, and plan to access their molecular cloning strategy to generate their vectors in order to obtain the sgRNA sequence for each of the genes targeted sufficient for Cell Ranger to assign perturbed genes to cells. Then, I will redo the Cell Ranger analysis using 10x's CIRSPR guide capture pipeline to assign perturbed genes to each cell.

- To replicate their figures by using their provided code, we downloaded the processed and normalized .h5ad single cell dataset to replicate the code in Perturbseq_GI (<https://github.com/thomasmaxwellnorman/Perturbseq_GI>) to hopefully replicate their analysis and final figures. We cloned their repository and created a virtual environment with the python version and libraries necessary to execute the code. Unfortunately, the .h5ad is quite large (65 gb), and trying to do a first-pass analysis with something like scanpy is too computationally expensive to perform on our MacBooks. We plan to generate a subsetted single-cell dataset using a subset of MTX gene-barcode counts provided from the authors. To process these MTX files (raw output from cell ranger) to call sgRNAs. We will use the scripts from the guide-calling (<https://github.com/josephreplogle/guide_calling>) codeset to call a sgRNA identity for each cell, then use the perturbseq_demo (<https://github.com/thomasmaxwellnorman/perturbseq_demo/tree/master>) codeset to generate and analyze a single-cell dataset.


Jander Kugelman: - Wrote transformers.R without tsne implementation -
About halfway through cell_population.R rewrite. This will likely be the
most problematic package to rewrite and test because of how long it is.

## Project Organization 
- figures we are attempting to recreate have
been moved to /images 
- progress on perturbseq package is in/src/perturbseqR/R

-   Cell Ranger output is stored in /processed/cell_ranger_output/20251107/
-   Scripts to recreate raw processing of data is stored in /src/subset_reanalysis
-   Scripts to recreate processing of normalized and filter single-cell dataset will be stored in /src/mtx_subset_analysis

## Struggles/Questions 
Jander: - I am still in the
process of rewriting the last script
(src/perturbseqR/R/cell_population.R). The only struggle I have right
now with it is that it is about 1000 lines of code, so it is taking me
awhile to make sure I am rewriting all of the functions correctly. - Dan
figured out that the tsne package for python used in the transformers.py
script in the perturbseq package is no longer functional, so I only
wrote the sklearn

Dan Le:
-   We unfortunately do not have the compute power on our MacBooks to effectively analyze the .h5ad single-cell dataset effectively. Would it be possible to gain access to the compute cluster, or would scaling back to an analysis of a subset (or an alternative single-cell dataset from the same publication but much smaller in size) be more feasible?

- Unfortunately, it seems like a lot of the metadata necessary to recreate the author's analyses have to be wrangled from what they actually provide. Additionally, some of the scripts the authors share on their repository have to be adjusted to be compatible with the datasets from the authors publication.


