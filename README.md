# InferCNV_Rscript
Place to put the script I used to process scRNA-seq CCRCC sample data before running inferCNV

## File needed to run thie scripts inlucdes
* h5 file for the scRNA-seq data (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices). Ensambl gene id and symbol mapping extracted from this file.
* Seurat Object after QC and ready to be analyzed saved as .rds file
* Genome position file. Can use the one in this repo "gencode_v19_gen_pos_sym_only.txt" or can be downloaded from TrinityCTAT
https://data.broadinstitute.org/Trinity/CTAT/cnv/
