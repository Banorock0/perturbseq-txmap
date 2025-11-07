#!/bin/bash

# This shell script contains the cell ranger CLI commands to upload the small subset of fastq files 
# from Repogle et al 2022


#txg files upload --project-id 599zYU7OiRDO84FGUzYCOyg /Users/cmdb/qb25-answers/perturbseq-txmap/rawdata/KD8_p1_0/p1_0_0_1_HLKJCDSXY

# upload to cellranger
for dir in p1_0_0_1_HLKJCDSXY p1_0_0_1_HLKKGDSXY p1_0_1_1_HLLCJDSXY p1_0_1_1_HLM2NDSXY; do txg --assumeyes files upload --project-id 599zYU7OiRDO84FGUzYCOyg /Users/cmdb/qb25-answers/perturbseq-txmap/rawdata/KD8_p1_0/$dir ; done 


