# Essential library
library(tidyverse)
library(cellrangerRkit)
library(Seurat)
library(infercnv)

# Setting up files pwd --------------------------------------------------------
folder_path   = 'TWAN-RCC/inferCNV_A1/'
ori_pos_path  = paste0(folder_path,'gencode_v19_gen_pos.complete.txt')
gen_pos_path  = paste0(folder_path,'gencode_v19_gen_pos_sym_only.txt')
annot_path    = paste0(folder_path,'A1_QC_annot.txt')
exp_path      = paste0(folder_path,'A1_QC_raw_expression.txt')
rds_path      = 'TWAN-RCC/ccRCC_A1_Aug2019.rds'
h5_path       = "TWAN-RCC/TWAN-RCC-D2/outs/raw_gene_bc_matrices_h5.h5" # h5 file for sample. Use to convert ensambl ID/Symbol
output_path   = paste0(folder_path,'inferCNV_result')

##################################################################
## Assisting function
## A function to toggle between Ensambl Id and symbol
input<- get_matrix_from_h5( h5_path )
ensl_id_sym <- fData(input)
input<- NULL

ensg_toggle <- function(input_marker){
  if(any(grepl('ENSG',input_marker))){output <- ensl_id_sym$symbol[match(input_marker,ensl_id_sym$id)]}
  else{output <- ensl_id_sym$id[match(input_marker,ensl_id_sym$symbol)]}
  output
}


##################################################################
# Process gene position file ----------------------------------------------
# Remake the file to have ensambl symbol only. Only need once if any.
tmp_pos <- read_tsv(ori_pos_path , col_name =F)
# # Keep ENSG ID
# name_tmp <- gsub('.*\\|',"",tmp_pos$X1)
# name_tmp2 <-gsub("\\..*",'',name_tmp)
# tmp_pos$X1 <- name_tmp2
# write.table(tmp_pos, 'gencode_v19_gen_pos_id_only.txt', col.names = F, quote =F, row.names = F, sep= '\t')
# tmp_pos <- read_tsv('gencode_v19_gen_pos_id_only.txt', col_name =F)

# Keep ENSG Symbol
tmp_pos$X1 <- gsub('\\|.*',"",tmp_pos$X1)
write.table(tmp_pos, gen_pos_path , col.names = F, quote =F, row.names = F, sep= '\t')

##################################################################
# Process from Seurat Object ----------------------------------------------
gene_pos_sym <- read.table(gen_pos_path ,sep='\t')
renal_a1 <- readRDS( rds_path )

# expression data - read
exp_df <- GetAssayData(renal_a1) %>% as.data.frame %>% rownames_to_column(var = 'gene')


tmp1<- ensg_toggle(exp_df$gene)
tmp_df <- data.frame(h5 = tmp1, biomart = tmp)
# expression data - sum duplicate genes
exp_df$gene <- ensg_toggle(exp_df$gene) # Transfer to ensnbol_sym
exp_df      <- exp_df %>% group_by(gene) %>% summarise_all(sum)

# expression data - remove gene not mapped in the gene position file
missing_sym  <- exp_df$gene[!exp_df$gene %in% gene_pos_sym$V1]
missing_row  <- match(missing_sym, exp_df$gene)
new_exp      <- exp_df[-missing_row,]
new_exp$gene %in% gene_pos_sym$V1 %>% all # Check if only gene with position left. Should be TRUE
new_exp$gene %>% duplicated %>% any # Check if duplicate still exist. Should be FALSE

# expression data - export
write.table(new_exp, file = exp_path, quote =F,sep='\t', row.names = F)

# Cell type annotation - export 
write.table(Idents(renal_a1) %>% as.data.frame, annot_path, quote= F, col.names = F, sep = '\t')

# Cell type annotation - test readback
annot_df <- read.table(annot_path, sep='\t')

# Output list of cell types
annot_df$V2 %>% unique %>% dput ## copy this and remove tumor related gene and put in 'nornal_groups'


##################################################################
# Main inferCNV function. Run if all files are ready ----------------------
# Assign normal_groups as reference
normal_groups <- c("NK(1)",'NK(2)','NK(3)','NK(4)'
                   ,"CD4Tcell",'CD8Tcell','Treg','Macrophage','DC',
                   'Mesangial','Bcell','Pro-Bcell','Endothelial','Astrocyte') # Replace this with normal cells in the sample

# Create inferCNV object --------------------------------------------------
infercnv_obj = CreateInfercnvObject(
  raw_counts_matrix= exp_path,
  annotations_file = annot_path,
  delim            = "\t",
  gene_order_file  = gen_pos_path,
  ref_group_names  = normal_groups)


# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff  = 0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir = output_path,  # dir is auto-created for storing outputs
                             cluster_by_groups = T,   # cluster
                             denoise = T,
                             HMM     = T
)


