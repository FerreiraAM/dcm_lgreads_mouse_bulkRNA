# rMATS data analysis

# Libraries
library("tidyverse")

# Functions
# Read files and add splice junction type
fun_read_SJ <- function(ls_file_names){
  lapply(ls_file_names, FUN = function(onefile){
    # Read file
    df <- read.table(onefile, header = TRUE)
    # Add splice junction type
    df$SJ_name <- unlist(strsplit(x = basename(onefile), "[.]"))[1]
    # Return
    df
  })
}
# Combine list of file in a data.frame
fun_comb_in_df <- function(ls_files, ls_file_names){
  # Rename
  names(ls_files) <- basename(ls_file_names)
  # Dataframe
  dplyr::bind_rows(ls_files)
}
# Basics stats from the files
fun_basic_stats <- function(df_files, fdr = 0.01){
  # Stats
  nb_event_per_SJtype <- df_files %>% 
    group_by(SJ_name) %>% 
    count(SJ_name)
  # Number of affected genes 
  nb_affected_genes <- length(unique(df_files$geneSymbol))
  # Number of gene affected by a significant SJ event
  df_files_sig <- dplyr::filter(df_files, FDR < fdr)
  nb_genes_sig <- length(unique(df_files_sig$geneSymbol))
  # Return
  list("nb_event_per_SJtype" = nb_event_per_SJtype,
       "nb_affected_genes" = nb_affected_genes,
       "nb_genes_sig" = nb_genes_sig)
}

# R636Q HOM vs WT
# List files
ls_R636Q_HOMvsWT_JC_filesname <- list.files(path = "rmats_directory/R636Q_HOMvsWT", 
                                pattern = "*JC.txt",
                                full.names = TRUE)
# Read files
ls_R636Q_HOMvsWT_JC_files <- fun_read_SJ(ls_R636Q_HOMvsWT_JC_filesname)
# Convert to a data.frame
df_R636Q_HOMvsWT_JC_files <- fun_comb_in_df(ls_R636Q_HOMvsWT_JC_files, 
                                            ls_R636Q_HOMvsWT_JC_filesname)
# Save file
write.table(df_R636Q_HOMvsWT_JC_files, 
            file = "rmats_directory/rmats_R636Q_HOMvsWT_JC_files.txt",
            col.names = TRUE,
            row.names = FALSE)
# Stats
fun_basic_stats(df_R636Q_HOMvsWT_JC_files)
# R636Q HET vs WT
# List files
ls_R636Q_HETvsWT_JC_filesname <- list.files(path = "rmats_directory/R636Q_HETvsWT", 
                                            pattern = "*JC.txt",
                                            full.names = TRUE)
# Read files
ls_R636Q_HETvsWT_JC_files <- fun_read_SJ(ls_R636Q_HETvsWT_JC_filesname)
# Convert to a data.frame
df_R636Q_HETvsWT_JC_files <- fun_comb_in_df(ls_R636Q_HETvsWT_JC_files, 
                                            ls_R636Q_HETvsWT_JC_filesname)
# Save file
write.table(df_R636Q_HETvsWT_JC_files, 
            file = "rmats_directory/rmats_R636Q_HETvsWT_JC_files.txt",
            col.names = TRUE,
            row.names = FALSE)
# Stats
fun_basic_stats(df_R636Q_HETvsWT_JC_files)

# P635L
# HOM vs WT
# List files
ls_P635L_HOMvsWT_JC_filesname <- list.files(path = "rmats_directory/P635L_HOMvsWT", 
                                            pattern = "*JC.txt",
                                            full.names = TRUE)
# Read files
ls_P635L_HOMvsWT_JC_files <- fun_read_SJ(ls_P635L_HOMvsWT_JC_filesname)
# Convert to a data.frame
df_P635L_HOMvsWT_JC_files <- fun_comb_in_df(ls_P635L_HOMvsWT_JC_files, 
                                            ls_P635L_HOMvsWT_JC_filesname)
# Save file
write.table(df_P635L_HOMvsWT_JC_files, 
            file = "rmats_directory/rmats_P635L_HOMvsWT_JC_files.txt",
            col.names = TRUE,
            row.names = FALSE)
# Stats
fun_basic_stats(df_P635L_HOMvsWT_JC_files)
# Overlap genes with DEXseq
df_P635L_JC_files_sig <- dplyr::filter(df_P635L_HOMvsWT_JC_files,
                                       FDR < 0.01)
## Load DEXseq results
P635L_dxd_res_wt_hom_gene_name_sig <- read.table("results/DEXSeq/table_P635L_dxd_res_wt_hom_gene_name_sig.txt",
                                                 header = TRUE)
v_P635L_dxd_res_wt_hom_gene_name_sig <- unique(P635L_dxd_res_wt_hom_gene_name_sig$gene_name)
## Intersect
inter_P635L_rmats_dexseq <- intersect(unique(df_P635L_JC_files_sig$geneSymbol), 
          v_P635L_dxd_res_wt_hom_gene_name_sig)
length(inter_P635L_rmats_dexseq)
