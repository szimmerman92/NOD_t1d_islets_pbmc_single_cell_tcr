Read me

NOD_t1d_islets_pbmc_single_cell_tcr is the GitHub repository that contains the code for the manuscript "Single-cell transcriptomics and TCR repertoire sequencing reveal the clonal expansion of islet infiltrating T cells in peripheral blood of NOD mice"

The published paper can be found here


The published datasets can be found at the following GEO repository: GSE200695


This directory has four subdirectories that can be from A to D:


- A.Sequence_processing
Initial sequence processing in cellranger pipeline

- B.Gene_expression_TCR_analysis
Runs scripts to match TCRs in both islets and PBMCs from NOD mice

- C.Pseudobulk_and_Trajectory_analysis
Runs pseubulk analysis for differential expression of genes between matching cells in blood and islets.
Run monocle3 analysis for trajectory inference.

- D.MachineLearning
Machine Learning analysis for Figure 4.