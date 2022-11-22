# rMATS data analysis

# Libraries
library("tidyverse")
library("pheatmap")

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

#### Load data ####
# Using the JC files that contains the junction counts

#### Initial experiment ####

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
# HET vs WT
# List files
ls_P635L_HETvsWT_JC_filesname <- list.files(path = "rmats_directory/P635L_HETvsWT", 
                                            pattern = "*JC.txt",
                                            full.names = TRUE)
# Read files
ls_P635L_HETvsWT_JC_files <- fun_read_SJ(ls_P635L_HETvsWT_JC_filesname)
# Convert to a data.frame
df_P635L_HETvsWT_JC_files <- fun_comb_in_df(ls_P635L_HETvsWT_JC_files, 
                                            ls_P635L_HETvsWT_JC_filesname)
# Save file
write.table(df_P635L_HETvsWT_JC_files, 
            file = "rmats_directory/rmats_P635L_HETvsWT_JC_files.txt",
            col.names = TRUE,
            row.names = FALSE)
# Stats
fun_basic_stats(df_P635L_HETvsWT_JC_files)
# # HET wo 3511 vs WT
# # List files
# ls_P635L_HETvsWT_wo3511_JC_filesname <- list.files(path = "rmats_directory/P635L_HETvsWT_wo3511", 
#                                                    pattern = "*JC.txt",
#                                                    full.names = TRUE)
# # Read files
# ls_P635L_HETvsWT_wo3511_JC_files <- fun_read_SJ(ls_P635L_HETvsWT_wo3511_JC_filesname)
# # Convert to a data.frame
# df_P635L_HETvsWT_wo3511_JC_files <- fun_comb_in_df(ls_P635L_HETvsWT_wo3511_JC_files, 
#                                                    ls_P635L_HETvsWT_wo3511_JC_filesname)
# # Save file
# write.table(df_P635L_HETvsWT_wo3511_JC_files, 
#             file = "rmats_directory/rmats_P635L_HET_wo3511vsWT_JC_files.txt",
#             col.names = TRUE,
#             row.names = FALSE)
# # Stats
# fun_basic_stats(df_P635L_HETvsWT_wo3511_JC_files)

#### Base editing experiments ####

# P635L Nterm_NRTH_Abe8e_and_Cterm_gRNA5 vs WT
# P635L_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT
# List files
ls_P635L_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_JC_filesname <- list.files(path = "rmats_directory/P635L_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT", 
                                                                                            pattern = "*JC.txt",
                                                                                            full.names = TRUE)
# Read files
ls_P635L_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_JC_files <- fun_read_SJ(ls_P635L_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_JC_filesname)
# Convert to a data.frame
df_P635L_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_files <- fun_comb_in_df(ls_P635L_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_JC_files, 
                                                                                         ls_P635L_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_JC_filesname)
# Save file
write.table(df_P635L_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_files, 
            file = "rmats_directory/rmats_P635L_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_JC_files.txt",
            col.names = TRUE,
            row.names = FALSE)
# Stats
fun_basic_stats(df_P635L_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_files)

# P635L_after_base_editing_Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT
# List files
ls_P635L_after_base_editing_Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_JC_filesname <- list.files(path = "rmats_directory/P635L_after_base_editing_Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT", 
                                                                                           pattern = "*JC.txt",
                                                                                           full.names = TRUE)
# Read files
ls_P635L_after_base_editing_Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_JC_files <- fun_read_SJ(ls_P635L_after_base_editing_Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_JC_filesname)
# Convert to a data.frame
df_P635L_after_base_editing_Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_files <- fun_comb_in_df(ls_P635L_after_base_editing_Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_JC_files, 
                                                                                        ls_P635L_after_base_editing_Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_JC_filesname)
# Save file
write.table(df_P635L_after_base_editing_Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_files, 
            file = "rmats_directory/rmats_P635L_after_base_editing_Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_JC_files.txt",
            col.names = TRUE,
            row.names = FALSE)
# Stats
fun_basic_stats(df_P635L_after_base_editing_Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_files)

# P635L_after_base_editing_PBSvsWT
# List files
ls_P635L_after_base_editing_PBSvsWT_JC_filesname <- list.files(path = "rmats_directory/P635L_after_base_editing_PBSvsWT", 
                                                               pattern = "*JC.txt",
                                                               full.names = TRUE)
# Read files
ls_P635L_after_base_editing_PBSvsWT_JC_files <- fun_read_SJ(ls_P635L_after_base_editing_PBSvsWT_JC_filesname)
# Convert to a data.frame
df_P635L_after_base_editing_PBSvsWT_files <- fun_comb_in_df(ls_P635L_after_base_editing_PBSvsWT_JC_files, 
                                                            ls_P635L_after_base_editing_PBSvsWT_JC_filesname)
# Save file
write.table(df_P635L_after_base_editing_PBSvsWT_files, 
            file = "rmats_directory/rmats_P635L_after_base_editing_PBSvsWT_JC_files.txt",
            col.names = TRUE,
            row.names = FALSE)
# Stats
fun_basic_stats(df_P635L_after_base_editing_PBSvsWT_files)

# R636Q_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6vsWT
# List files
ls_R636Q_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6vsWT_JC_filesname <- list.files(path = "rmats_directory/R636Q_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6vsWT", 
                                                                                           pattern = "*JC.txt",
                                                                                           full.names = TRUE)
# Read files
ls_R636Q_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6vsWT_JC_files <- fun_read_SJ(ls_R636Q_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6vsWT_JC_filesname)
# Convert to a data.frame
df_R636Q_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6vsWT_files <- fun_comb_in_df(ls_R636Q_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6vsWT_JC_files, 
                                                                                        ls_R636Q_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6vsWT_JC_filesname)
# Save file
write.table(df_R636Q_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6vsWT_files, 
            file = "rmats_directory/rmats_R636Q_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6vsWT_JC_files.txt",
            col.names = TRUE,
            row.names = FALSE)
# Stats
fun_basic_stats(df_R636Q_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6vsWT_files)

# R636Q_after_base_editing_PBSvsWT
# List files
ls_R636Q_after_base_editing_PBSvsWT_JC_filesname <- list.files(path = "rmats_directory/R636Q_after_base_editing_PBSvsWT", 
                                                               pattern = "*JC.txt",
                                                               full.names = TRUE)
# Read files
ls_R636Q_after_base_editing_PBSvsWT_JC_files <- fun_read_SJ(ls_R636Q_after_base_editing_PBSvsWT_JC_filesname)
# Convert to a data.frame
df_R636Q_after_base_editing_PBSvsWT_files <- fun_comb_in_df(ls_R636Q_after_base_editing_PBSvsWT_JC_files, 
                                                            ls_R636Q_after_base_editing_PBSvsWT_JC_filesname)
# Save file
write.table(df_R636Q_after_base_editing_PBSvsWT_files, 
            file = "rmats_directory/rmats_R636Q_after_base_editing_PBSvsWT_JC_files.txt",
            col.names = TRUE,
            row.names = FALSE)
# Stats
fun_basic_stats(df_R636Q_after_base_editing_PBSvsWT_files)

# R636Q_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6vsPBS
# List files
ls_R636Q_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6vsPBS_JC_filesname <- list.files(path = "rmats_directory/R636Q_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6vsPBS", 
                                                               pattern = "*JC.txt",
                                                               full.names = TRUE)
# Read files
ls_R636Q_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6vsPBS_JC_files <- fun_read_SJ(ls_R636Q_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6vsPBS_JC_filesname)
# Convert to a data.frame
df_R636Q_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6vsPBS_files <- fun_comb_in_df(ls_R636Q_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6vsPBS_JC_files, 
                                                            ls_R636Q_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6vsPBS_JC_filesname)
# Save file
write.table(df_R636Q_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6vsPBS_files, 
            file = "rmats_directory/rmats_R636Q_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6vsPBS_JC_files.txt",
            col.names = TRUE,
            row.names = FALSE)
# Stats
fun_basic_stats(df_R636Q_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6vsPBS_files)

# P635L_after_base_editing_PBSvsNterm_SpRY_and_Cterm_SpRY_gRNA5
# List files
ls_P635L_after_base_editing_PBSvsNterm_SpRY_and_Cterm_SpRY_gRNA5_JC_filesname <- list.files(path = "rmats_directory/P635L_after_base_editing_PBSvsNterm_SpRY_and_Cterm_SpRY_gRNA5", 
                                                                                                  pattern = "*JC.txt",
                                                                                                  full.names = TRUE)
# Read files
ls_P635L_after_base_editing_PBSvsNterm_SpRY_and_Cterm_SpRY_gRNA5_JC_files <- fun_read_SJ(ls_P635L_after_base_editing_PBSvsNterm_SpRY_and_Cterm_SpRY_gRNA5_JC_filesname)
# Convert to a data.frame
df_P635L_after_base_editing_PBSvsNterm_SpRY_and_Cterm_SpRY_gRNA5_files <- fun_comb_in_df(ls_P635L_after_base_editing_PBSvsNterm_SpRY_and_Cterm_SpRY_gRNA5_JC_files, 
                                                                                               ls_P635L_after_base_editing_PBSvsNterm_SpRY_and_Cterm_SpRY_gRNA5_JC_filesname)
# Save file
write.table(df_P635L_after_base_editing_PBSvsNterm_SpRY_and_Cterm_SpRY_gRNA5_files, 
            file = "rmats_directory/rmats_P635L_after_base_editing_PBSvsNterm_SpRY_and_Cterm_SpRY_gRNA5_JC_files.txt",
            col.names = TRUE,
            row.names = FALSE)
# Stats
fun_basic_stats(df_P635L_after_base_editing_PBSvsNterm_SpRY_and_Cterm_SpRY_gRNA5_files)

# P635L_after_base_editing_PBSvsNterm_NRTH_Abe8e_and_Cterm_gRNA5
# List files
ls_P635L_after_base_editing_PBSvsNterm_NRTH_Abe8e_and_Cterm_gRNA5_JC_filesname <- list.files(path = "rmats_directory/P635L_after_base_editing_PBSvsNterm_NRTH_Abe8e_and_Cterm_gRNA5", 
                                                                                            pattern = "*JC.txt",
                                                                                            full.names = TRUE)
# Read files
ls_P635L_after_base_editing_PBSvsNterm_NRTH_Abe8e_and_Cterm_gRNA5_JC_files <- fun_read_SJ(ls_P635L_after_base_editing_PBSvsNterm_NRTH_Abe8e_and_Cterm_gRNA5_JC_filesname)
# Convert to a data.frame
df_P635L_after_base_editing_PBSvsNterm_NRTH_Abe8e_and_Cterm_gRNA5_files <- fun_comb_in_df(ls_P635L_after_base_editing_PBSvsNterm_NRTH_Abe8e_and_Cterm_gRNA5_JC_files, 
                                                                                         ls_P635L_after_base_editing_PBSvsNterm_NRTH_Abe8e_and_Cterm_gRNA5_JC_filesname)
# Save file
write.table(df_P635L_after_base_editing_PBSvsNterm_NRTH_Abe8e_and_Cterm_gRNA5_files, 
            file = "rmats_directory/rmats_P635L_after_base_editing_PBSvsNterm_NRTH_Abe8e_and_Cterm_gRNA5_JC_files.txt",
            col.names = TRUE,
            row.names = FALSE)
# Stats
fun_basic_stats(df_P635L_after_base_editing_PBSvsNterm_NRTH_Abe8e_and_Cterm_gRNA5_files)

# P635L_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsNterm_SpRY_and_Cterm_SpRY_gRNA5
# List files
ls_P635L_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsNterm_SpRY_and_Cterm_SpRY_gRNA5_JC_filesname <- list.files(path = "rmats_directory/P635L_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsNterm_SpRY_and_Cterm_SpRY_gRNA5", 
                                                                                             pattern = "*JC.txt",
                                                                                             full.names = TRUE)
# Read files
ls_P635L_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsNterm_SpRY_and_Cterm_SpRY_gRNA5_JC_files <- fun_read_SJ(ls_P635L_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsNterm_SpRY_and_Cterm_SpRY_gRNA5_JC_filesname)
# Convert to a data.frame
df_P635L_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsNterm_SpRY_and_Cterm_SpRY_gRNA5_files <- fun_comb_in_df(ls_P635L_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsNterm_SpRY_and_Cterm_SpRY_gRNA5_JC_files, 
                                                                                          ls_P635L_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsNterm_SpRY_and_Cterm_SpRY_gRNA5_JC_filesname)
# Save file
write.table(df_P635L_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsNterm_SpRY_and_Cterm_SpRY_gRNA5_files, 
            file = "rmats_directory/rmats_P635L_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsNterm_SpRY_and_Cterm_SpRY_gRNA5_JC_files.txt",
            col.names = TRUE,
            row.names = FALSE)
# Stats
fun_basic_stats(df_P635L_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsNterm_SpRY_and_Cterm_SpRY_gRNA5_files)

###### Overlap ########
# Previously --> replace by tidyverse command below
# # Run in the script rmats_data_analysis_cluster.R
# ls_sameevent_pergene <- readRDS(file = "rmats_directory/ls_sameevent_pergene.rds")
# df_sameevent_pergene <- bind_rows(ls_sameevent_pergene)
# df_sameevent_pergene_sig <- df_sameevent_pergene %>% 
#   group_by(geneSymbol, SJ_name, SJ_ID) %>% 
#   dplyr::filter((all(FDR[comparison == "HOMvsWT"] < 0.01) | all(FDR[comparison == "HOMvsWT"] == 0)) &
#                   (all(IncLevelDifference[comparison == "HOMvsWT"] < -0.1) | 
#                      all(IncLevelDifference[comparison == "HOMvsWT"] > 0.1)))
# HOM vs WT
# Add a column in each dataframe
df_P635L_HOMvsWT_JC_files$mutant <- "P635L"
df_P635L_HOMvsWT_JC_files$comparison <- "HOMvsWT"
df_R636Q_HOMvsWT_JC_files$mutant <- "R636Q"
df_R636Q_HOMvsWT_JC_files$comparison <- "HOMvsWT"
df_P635L_HETvsWT_JC_files$mutant <- "P635L"
df_P635L_HETvsWT_JC_files$comparison <- "HETvsWT"
df_R636Q_HETvsWT_JC_files$mutant <- "R636Q"
df_R636Q_HETvsWT_JC_files$comparison <- "HETvsWT"
# Combine
df_all_JC_files <- rbind(df_P635L_HOMvsWT_JC_files,
                         df_R636Q_HOMvsWT_JC_files,
                         df_P635L_HETvsWT_JC_files,
                         df_R636Q_HETvsWT_JC_files)
# Add SJ ID
df_all_JC_files <- df_all_JC_files %>% 
  group_by(geneSymbol, SJ_name, longExonStart_0base, longExonEnd,
           shortES, shortEE, flankingES, flankingEE,
           X1stExonStart_0base, X1stExonEnd, X2ndExonStart_0base, X2ndExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE,
           riExonEnd, exonStart_0base, exonEnd) %>% 
  mutate(SJ_ID = cur_group_id())
# Filter for SJ events found more in all 4 
df_all_JC_files_SJinall4 <- df_all_JC_files %>% 
  group_by(SJ_ID) %>% 
  dplyr::filter(n() == 4)
# Filter the significant genes based on HOM vs WT comparison
# Event with FDR < 0.01 or = 0 
# FDR and p-value = 0 
# https://groups.google.com/g/rmats-user-group/c/TW534af62fg
# Zero P-value and FDR are generated when the actual value is smaller than the numerical accuracy cutoff, 
# which is usually 2.2e-16. Therefore, zero P-values can be interpreted as P<= 2.2e-16.
# And IncLevelDifference < -0.1 or > 0.1 
df_sameevent_pergene_sig <- df_all_JC_files_SJinall4 %>% 
  group_by(SJ_ID) %>% 
  dplyr::filter((all(FDR[comparison == "HOMvsWT"] < 0.01) | all(FDR[comparison == "HOMvsWT"] == 0)) &
                  (all(IncLevelDifference[comparison == "HOMvsWT"] < -0.1) | 
                     all(IncLevelDifference[comparison == "HOMvsWT"] > 0.1)))
# Number of genes
length(unique(df_sameevent_pergene_sig$geneSymbol))
# Write the results
write.table(df_sameevent_pergene_sig,
            file = "rmats_directory/table_allSJs_common_siginHOMvsWT.txt",
            row.names = FALSE)

# Heatmap
# All significant SJ common in HOM vs WT comparisons found in the 4 comparisons
# Format
# Filter values
df_sameevent_pergene_sig_filtered <- dplyr::select(df_sameevent_pergene_sig,
                                                   geneSymbol, SJ_name, SJ_ID, mutant, comparison, IncLevelDifference)
# Combine columns
df_sameevent_pergene_sig_filtered_tomat <- unite(df_sameevent_pergene_sig_filtered, 
                                                 gene_SJ, 
                                                 geneSymbol, SJ_name, SJ_ID)
df_sameevent_pergene_sig_filtered_tomat <- unite(df_sameevent_pergene_sig_filtered_tomat, 
                                                 mutant_comparison, 
                                                 mutant, comparison)
# Pivot wider
df_sameevent_pergene_sig_filtered_tomat_w <- pivot_wider(df_sameevent_pergene_sig_filtered_tomat,
                                                         names_from = mutant_comparison,
                                                         values_from = IncLevelDifference)
# Rownames
mat_sameevent_pergene_sig_filtered_w <- column_to_rownames(df_sameevent_pergene_sig_filtered_tomat_w,
                                                           var = "gene_SJ")
# Plot
pdf("rmats_directory/heatmap_rmats_allSJs_common_siginHOMvsWT.pdf", width = 6, height = 7)
pheatmap(mat_sameevent_pergene_sig_filtered_w,
         fontsize = 5,
         main = "All significant SJs common in HOM vs WT comparisons found in the 4 comparisons")
dev.off()

# P635L HET, P635L HOM, R636Q HET, R636Q HOM
pdf("rmats_directory/DRAFT_heatmap_rmats_allSJs_common_siginHOMvsWT.pdf", width = 11, height = 3)
pheatmap(t(mat_sameevent_pergene_sig_filtered_w[,c("P635L_HETvsWT",
                                                   "P635L_HOMvsWT",
                                                   "R636Q_HETvsWT",
                                                   "R636Q_HOMvsWT")]),
         color=colorRampPalette(c("navy", "white", "red"))(50),
         fontsize_col = 7,
         cluster_rows = FALSE,
         labels_row = c("P635L HET", "P635L HOM", "R636Q HET", "R636Q HOM"),
         labels_col = sub("_.*", "", rownames(mat_sameevent_pergene_sig_filtered_w)))
dev.off()

# Add SJ not overlaping but important genes
# Cryab - Hdnr - Hectd2os - Mia2 - Ryr2 - Slmap - Usp48
# Intersect all
intersect_all_genesymbol <- Reduce(intersect, list(unique(df_P635L_HOMvsWT_JC_files$geneSymbol),
                                                   unique(df_R636Q_HOMvsWT_JC_files$geneSymbol),
                                                   unique(df_P635L_HETvsWT_JC_files$geneSymbol),
                                                   unique(df_R636Q_HETvsWT_JC_files$geneSymbol)))
# Filtering
# Common gene symbols
df_all_JC_files_common_genesymbols <- dplyr::filter(df_all_JC_files, 
                                                    geneSymbol %in% intersect_all_genesymbol)
# Filter the one missing
miss_genes <- c("Cryab", "Hdnr", "Hectd2os", "Mia2", "Ryr2", "Usp48")
df_miss_genes_sig <- dplyr::filter(df_all_JC_files_common_genesymbols, geneSymbol %in% miss_genes) %>% 
  group_by(geneSymbol, SJ_name, longExonStart_0base, longExonEnd,
           shortES, shortEE, flankingES, flankingEE,
           X1stExonStart_0base, X1stExonEnd, X2ndExonStart_0base, X2ndExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE,
           riExonEnd, exonStart_0base, exonEnd) %>% 
  dplyr::filter((all(FDR[comparison == "HOMvsWT"] < 0.01) | all(FDR[comparison == "HOMvsWT"] == 0)) &
                  (all(IncLevelDifference[comparison == "HOMvsWT"] < -0.1) | 
                     all(IncLevelDifference[comparison == "HOMvsWT"] > 0.1)))
# Add SJ IDs
df_miss_genes_sig <- df_miss_genes_sig %>% 
  group_by(geneSymbol, SJ_name, longExonStart_0base, longExonEnd,
           shortES, shortEE, flankingES, flankingEE,
           X1stExonStart_0base, X1stExonEnd, X2ndExonStart_0base, X2ndExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE,
           riExonEnd, exonStart_0base, exonEnd) %>% 
  mutate(SJ_ID = cur_group_id())
# Filter for SJ events found more than once
df_miss_genes_sig <- df_miss_genes_sig %>% 
  group_by(SJ_ID) %>% 
  dplyr::filter(n() > 1)
# Format
# Filter values
df_miss_genes_sig_filtered <- df_miss_genes_sig[,c("geneSymbol", "SJ_name", "SJ_ID", "mutant", 
                                                  "comparison", "IncLevelDifference")]
# Combine columns
df_miss_genes_sig_filtered_tomat <- unite(df_miss_genes_sig_filtered, 
                                                 gene_SJ, 
                                                 geneSymbol, SJ_name, SJ_ID)
df_miss_genes_sig_filtered_tomat <- unite(df_miss_genes_sig_filtered_tomat, 
                                                 mutant_comparison, 
                                                 mutant, comparison)
# Pivot wider
df_miss_genes_sig_filtered_tomat_w <- pivot_wider(df_miss_genes_sig_filtered_tomat,
                                                         names_from = mutant_comparison,
                                                         values_from = IncLevelDifference)
# Rownames
mat_miss_genes_sig_sig_filtered_w <- column_to_rownames(df_miss_genes_sig_filtered_tomat_w,
                                                           var = "gene_SJ")
# Combine data
comb_SJevent_withmissgenes_sig_filtered_w <- bind_rows(mat_sameevent_pergene_sig_filtered_w,
      mat_miss_genes_sig_sig_filtered_w)
pdf("rmats_directory/VF_heatmap_rmats_allSJs_common_and_missinggenes_siginHOMvsWT.pdf", 
    width = 4, height = 11)
pheatmap(comb_SJevent_withmissgenes_sig_filtered_w[,c("P635L_HETvsWT",
                                                   "P635L_HOMvsWT",
                                                   "R636Q_HETvsWT",
                                                   "R636Q_HOMvsWT")],
         color=colorRampPalette(c("navy", "white", "red"))(50),
         fontsize = 13,
         cluster_cols = FALSE,
         labels_col = c("P635L HET", "P635L HOM", "R636Q HET", "R636Q HOM"),
         labels_row = sub("_.*", "", rownames(comb_SJevent_withmissgenes_sig_filtered_w)),
         na_col = "gray58")
dev.off()

#### After base editing samples ####

#### R636Q ####
# identifying the significant SJ events between PBS vs WT for R636Q 
# PBSvsWT - Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT
# Add a column in each dataframe
df_R636Q_after_base_editing_PBSvsWT_files$mutant <- "R636Q_after_base_editing"
df_R636Q_after_base_editing_PBSvsWT_files$comparison <- "PBSvsWT"
df_R636Q_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6vsWT_files$mutant <- "R636Q_after_base_editing"
df_R636Q_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6vsWT_files$comparison <- "Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6"
# Combine
df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6 <- rbind(df_R636Q_after_base_editing_PBSvsWT_files, 
                                                                              df_R636Q_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6vsWT_files)
# Add SJ ID
df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6 <- df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6 %>% 
  group_by(geneSymbol, SJ_name, longExonStart_0base, longExonEnd,
           shortES, shortEE, flankingES, flankingEE,
           X1stExonStart_0base, X1stExonEnd, X2ndExonStart_0base, X2ndExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE,
           riExonEnd, exonStart_0base, exonEnd) %>% 
  mutate(SJ_ID = cur_group_id())
# Filter for SJ events found more than once
df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1 <- df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6 %>% 
  group_by(SJ_ID) %>% 
  dplyr::filter(n() > 1)
# Significant
df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_sig <- df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1 %>% 
  group_by(SJ_ID) %>% 
  dplyr::filter((all(FDR[comparison == "PBSvsWT"] < 0.01) | all(FDR[comparison == "PBSvsWT"] == 0)) &
                  (all(IncLevelDifference[comparison == "PBSvsWT"] < -0.1) | 
                     all(IncLevelDifference[comparison == "PBSvsWT"] > 0.1)))
# Order factor
df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_sig <- df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_sig %>% 
  mutate(comparison = factor(comparison, levels = c("PBSvsWT", "Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6")))
# Write
write.table(df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_sig,
            file = "rmats_directory/table_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_sig.txt",
            row.names = FALSE)
# Wide format for Markus
df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_sig_wd <- df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_sig %>% 
  dplyr::select(geneSymbol, SJ_ID, IncLevelDifference, comparison) %>% 
  pivot_wider(names_from = comparison, values_from = IncLevelDifference)
write.table(df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_sig_wd,
            file = "rmats_directory/table_df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_sig_wd.txt",
            row.names = FALSE)
# SJ events classification
df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_sig_wd_delta <- 
  df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_sig_wd %>% 
  mutate(diff = PBSvsWT - Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6)
# # Classification function --> does not work; dplyr BUGGGGGG
# fun_classification_SJevents <- function(df_data, x_original, x_edited){
#   df_class <- df_data %>% 
#     dplyr::rowwise() %>% 
#     # mutate(SJ_type = if_else((x_edited) <= 0, "rescue", "nope"))
#   dplyr::mutate(SJ_type = case_when(abs(diff) > 0.1 & (x_original) > 0 & (x_edited) >= 0 & (x_original) > (x_edited) ~ "rescue",
#                            abs(diff) > 0.1 & (x_original) < 0 & (x_edited) <= 0 & (x_original) < (x_edited) ~ "rescue",
#                            abs(diff) > 0.1 & (x_original) > 0 & (x_edited) >= 0 & (x_original) < (x_edited) ~ "reversal",
#                            abs(diff) > 0.1 & (x_original) < 0 & (x_edited) <= 0 & (x_original) > (x_edited) ~ "reversal",
#                            abs(diff) > 0.1 & (x_original) > 0 & (x_original) > (x_edited) &
#                              ((x_edited) >= -0.2 & (x_edited) <= 0.2) ~ "rescue",
#                            abs(diff) > 0.1 & (x_original) > 0 & (x_edited) < -0.2 ~ "reversal",
#                            abs(diff) > 0.1 & (x_original) < 0 & (x_original) < (x_edited) &
#                              ((x_edited) >= -0.2 & (x_edited) <= 0.2) ~ "rescue",
#                            abs(diff) > 0.1 & (x_original) < 0 & (x_edited) > 0.2) ~ "reversal"))
#   # Replace NA
#   replace_na(df_class, list(SJ_type = "diff_inf_0.1"))
# }
# R636Q
df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_sig_wd_delta <- 
  df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_sig_wd_delta %>% 
  mutate(SJ_type = case_when(abs(diff) > 0.1 & PBSvsWT > 0 & Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6 >= 0 & PBSvsWT > Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6 ~ "rescue",
                             abs(diff) > 0.1 & PBSvsWT < 0 & Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6 <= 0 & PBSvsWT < Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6 ~ "rescue",
                             abs(diff) > 0.1 & PBSvsWT > 0 & Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6 >= 0 & PBSvsWT < Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6 ~ "reversal",
                             abs(diff) > 0.1 & PBSvsWT < 0 & Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6 <= 0 & PBSvsWT > Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6 ~ "reversal",
                             abs(diff) > 0.1 & PBSvsWT > 0 & 
                               PBSvsWT > (Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6) & 
                               (Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6 >= -0.2 & Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6 <= 0.2) ~ "rescue",
                             abs(diff) > 0.1 & PBSvsWT > 0 & 
                               # PBSvsWT > (Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6) & 
                               # (Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6 < -0.2 | Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6 > 0.2) ~ "reversal",
                               Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6 < -0.2 ~ "reversal",
                             abs(diff) > 0.1 & PBSvsWT < 0 & 
                               (PBSvsWT) < Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6 & 
                               (Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6 >= -0.2 & Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6 <= 0.2) ~ "rescue",
                             abs(diff) > 0.1 & PBSvsWT < 0 &  
                               # (PBSvsWT) < Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6 & 
                               # (Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6 < -0.2 | Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6 > 0.2) ~ "reversal"))
                               Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6 > 0.2 ~ "reversal"))

# Replace NA 
df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_sig_wd_delta <- 
  df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_sig_wd_delta %>% 
  replace_na(list(SJ_type = "diff_inf_0.1"))
# long format to plot
df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_sig_lg_delta <- pivot_longer(df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_sig_wd_delta,
             cols = c(PBSvsWT, Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6),
             names_to = "comparison",
             values_to = "IncLevelDifference")
df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_sig_lg_delta <- 
  df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_sig_lg_delta %>% 
  mutate(comparison = factor(comparison, levels = c("PBSvsWT", "Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6")))
# Add Ttn column
df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_sig_lg_delta_ttn <- 
  df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_sig_lg_delta %>% 
  mutate("Ttn" = geneSymbol == "Ttn")
# Add color variable
df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_sig_lg_delta_ttn <- 
  df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_sig_lg_delta_ttn %>% 
  mutate(plot_color = if_else(Ttn == TRUE, true = "Ttn", false = as.character(SJ_type)))
# Factor order
df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_sig_lg_delta_ttn <- 
  df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_sig_lg_delta_ttn %>% 
  mutate(SJ_type = factor(SJ_type, levels = c("rescue", "reversal", "diff_inf_0.1")))
# Plot
pdf("rmats_directory/VF_R636Q_after_base_editing_dotlineplot_allSJs_common_sigin_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_deltadeltaPSI_facet.pdf", 
    width = 4, height = 6)
ggplot(df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_sig_lg_delta_ttn,
       aes(x = comparison, y = IncLevelDifference, group = SJ_ID, color = plot_color)) +
  geom_point(alpha = 0.6, aes(size = plot_color)) +
  geom_line(alpha = 0.6, aes(size = plot_color)) + 
  facet_wrap(vars(SJ_type), strip.position="right", ncol = 1, labeller = as_labeller(c("rescue" = "Rescue",
                                                             "reversal" = "Mis-spliced", 
                                                             "diff_inf_0.1" = "Unchanged"))) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  scale_color_manual(name = NULL, 
                     values = c("diff_inf_0.1" = "gray58", "rescue" = "red", "reversal" = "gray30", "Ttn" = "blue"),
                     breaks = c("Ttn")) +
  scale_size_manual(guide = 'none', values = c(0.2, 0.2, 0.2, 0.8)) +
  ggtitle("R636Q") +
  ylab("Delta delta PSI") +
  scale_x_discrete(name = element_blank(), 
                   labels = c("PBS", "Abe8e")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  annotate("rect", xmin = 0, xmax = 3, ymin = -0.2, ymax = 0.2, 
           alpha = .2) +
  guides(colour = guide_legend(override.aes = list(alpha = 1,
                                                   linetype = NULL,
                                                   size = 2))) 
dev.off()

#### P635L ####
# PBSvsWT - Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT
# Add a column in each dataframe
df_P635L_after_base_editing_PBSvsWT_files$mutant <- "P635L_after_base_editing"
df_P635L_after_base_editing_PBSvsWT_files$comparison <- "PBSvsWT"
df_P635L_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_files$mutant <- "P635L_after_base_editing"
df_P635L_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_files$comparison <- "Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT"
# Combine
df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5 <- rbind(df_P635L_after_base_editing_PBSvsWT_files, 
                                                                              df_P635L_after_base_editing_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_files)
# Add SJ ID
df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5 <- df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5 %>% 
  group_by(geneSymbol, SJ_name, longExonStart_0base, longExonEnd,
           shortES, shortEE, flankingES, flankingEE,
           X1stExonStart_0base, X1stExonEnd, X2ndExonStart_0base, X2ndExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE,
           riExonEnd, exonStart_0base, exonEnd) %>% 
  mutate(SJ_ID = cur_group_id())
# Filter for SJ events found more than once
df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1 <- df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5 %>% 
  group_by(SJ_ID) %>% 
  dplyr::filter(n() > 1)
# Significant
df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_sig <- df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1 %>% 
  group_by(SJ_ID) %>% 
  dplyr::filter((all(FDR[comparison == "PBSvsWT"] < 0.01) | all(FDR[comparison == "PBSvsWT"] == 0)) &
                  (all(IncLevelDifference[comparison == "PBSvsWT"] < -0.1) | 
                     all(IncLevelDifference[comparison == "PBSvsWT"] > 0.1)))
# Order factor
df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_sig <- df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_sig %>% 
  mutate(comparison = factor(comparison, levels = c("PBSvsWT", "Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT")))
# Write
write.table(df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_sig,
            file = "rmats_directory/table_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_sig.txt",
            row.names = FALSE)
# Wide format for Markus
df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_sig_wd <- df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_sig %>% 
  dplyr::select(geneSymbol, SJ_ID, IncLevelDifference, comparison) %>% 
  pivot_wider(names_from = comparison, values_from = IncLevelDifference)
write.table(df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_sig_wd,
            file = "rmats_directory/table_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_sig_wd.txt",
            row.names = FALSE)
# SJ events classification
df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_sig_wd_delta <- 
  df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_sig_wd %>% 
  mutate(diff = PBSvsWT - Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT) 
# P635L
df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_sig_wd_delta <- 
  df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_sig_wd_delta %>% 
  mutate(SJ_type = case_when(abs(diff) > 0.1 & PBSvsWT > 0 & Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT >= 0 & PBSvsWT > Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT ~ "rescue",
                             abs(diff) > 0.1 & PBSvsWT < 0 & Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT <= 0 & PBSvsWT < Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT ~ "rescue",
                             abs(diff) > 0.1 & PBSvsWT > 0 & Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT >= 0 & PBSvsWT < Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT ~ "reversal",
                             abs(diff) > 0.1 & PBSvsWT < 0 & Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT <= 0 & PBSvsWT > Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT ~ "reversal",
                             abs(diff) > 0.1 & PBSvsWT > 0 &  
                               PBSvsWT > (Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT) & 
                               (Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT >= -0.2 & Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT <= 0.2) ~ "rescue",
                             abs(diff) > 0.1 & PBSvsWT > 0 & 
                               PBSvsWT > (Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT) & 
                               (Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT < -0.2 | Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT > 0.2) ~ "reversal",
                             abs(diff) > 0.1 & PBSvsWT < 0 & 
                               (PBSvsWT) < Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT & 
                               (Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT >= -0.2 & Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT <= 0.2) ~ "rescue",
                             abs(diff) > 0.1 & PBSvsWT < 0 & 
                               (PBSvsWT) < Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT & 
                               (Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT < -0.2 | Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT > 0.2) ~ "reversal"))
# Replace NA 
df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_sig_wd_delta <- 
  df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_sig_wd_delta %>% 
  replace_na(list(SJ_type = "diff_inf_0.1"))
# long format to plot
df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_SJmore1_sig_lg_delta <- pivot_longer(df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_sig_wd_delta,
                                                                                                               cols = c(PBSvsWT, Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT),
                                                                                                               names_to = "comparison",
                                                                                                               values_to = "IncLevelDifference")
df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_SJmore1_sig_lg_delta <- 
  df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_SJmore1_sig_lg_delta %>% 
  mutate(comparison = factor(comparison, levels = c("PBSvsWT", "Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT")))
# Add Ttn column
df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_SJmore1_sig_lg_delta_ttn <- 
  df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_SJmore1_sig_lg_delta %>% 
  mutate("Ttn" = geneSymbol == "Ttn")
# Add color variable
df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_SJmore1_sig_lg_delta_ttn <- 
  df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_SJmore1_sig_lg_delta_ttn %>% 
  mutate(plot_color = if_else(Ttn == TRUE, true = "Ttn", false = as.character(SJ_type)))
# Factor order
df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_SJmore1_sig_lg_delta_ttn <- 
  df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_SJmore1_sig_lg_delta_ttn %>% 
  mutate(SJ_type = factor(SJ_type, levels = c("rescue", "reversal", "diff_inf_0.1")))
# Plot
pdf("rmats_directory/VF_P635L_after_base_editing_dotlineplot_allSJs_common_sigin_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_deltadeltaPSI_facet.pdf", 
    width = 4, height = 6)
ggplot(df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_SJmore1_sig_lg_delta_ttn,
       aes(x = comparison, y = IncLevelDifference, group = SJ_ID, color = plot_color)) +
  geom_point(alpha = 0.6, aes(size = plot_color)) +
  geom_line(alpha = 0.6, aes(size = plot_color)) + 
  facet_wrap(vars(SJ_type), strip.position="right", ncol = 1, labeller = as_labeller(c("rescue" = "Rescue",
                                                                                       "reversal" = "Mis-spliced", 
                                                                                       "diff_inf_0.1" = "Unchanged"))) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  scale_color_manual(name = NULL, 
                     values = c("diff_inf_0.1" = "gray58", "rescue" = "red", "reversal" = "gray30", "Ttn" = "blue"),
                     breaks = c("Ttn")) +
  scale_size_manual(guide = 'none', values = c(0.2, 0.2, 0.2, 0.8)) +
  ggtitle("P635L") +
  ylab("Delta delta PSI") +
  scale_x_discrete(name = element_blank(), 
                   labels = c("PBS", "Abe8e")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  annotate("rect", xmin = 0, xmax = 3, ymin = -0.2, ymax = 0.2, 
           alpha = .2) +
  guides(colour = guide_legend(override.aes = list(alpha = 1,
                                                   linetype = NULL,
                                                   size = 2))) 
dev.off()

# PBSvsWT - Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT
# Add column
df_P635L_after_base_editing_Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_files$mutant <- "P635L_after_base_editing"
df_P635L_after_base_editing_Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_files$comparison <- "Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT"
# Combine
df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5 <- rbind(df_P635L_after_base_editing_PBSvsWT_files, 
                                                                             df_P635L_after_base_editing_Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_files)
# Add SJ ID
df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5 <- df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5 %>% 
  group_by(geneSymbol, SJ_name, longExonStart_0base, longExonEnd,
           shortES, shortEE, flankingES, flankingEE,
           X1stExonStart_0base, X1stExonEnd, X2ndExonStart_0base, X2ndExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE,
           riExonEnd, exonStart_0base, exonEnd) %>% 
  mutate(SJ_ID = cur_group_id())
# Filter for SJ events found more than once
df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1 <- df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5 %>% 
  group_by(SJ_ID) %>% 
  dplyr::filter(n() > 1)
# Significant
df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_sig <- df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1 %>% 
  group_by(SJ_ID) %>% 
  dplyr::filter((all(FDR[comparison == "PBSvsWT"] < 0.01) | all(FDR[comparison == "PBSvsWT"] == 0)) &
                  (all(IncLevelDifference[comparison == "PBSvsWT"] < -0.1) | 
                     all(IncLevelDifference[comparison == "PBSvsWT"] > 0.1)))
# Order factor
df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_sig <- df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_sig %>% 
  mutate(comparison = factor(comparison, levels = c("PBSvsWT", "Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT")))
# Write
write.table(df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_sig,
            file = "rmats_directory/table_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_sig.txt",
            row.names = FALSE)
# Wide format for Markus
df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_sig_wd <- df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_sig %>% 
  dplyr::select(geneSymbol, SJ_ID, IncLevelDifference, comparison) %>% 
  pivot_wider(names_from = comparison, values_from = IncLevelDifference)
write.table(df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_sig_wd,
            file = "rmats_directory/table_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_sig_wd.txt",
            row.names = FALSE)
# SJ events classification
df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_sig_wd_delta <- 
  df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_sig_wd %>% 
  mutate(diff = PBSvsWT - Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT) 
# P635L
df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_sig_wd_delta <- 
  df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_sig_wd_delta %>% 
  mutate(SJ_type = case_when(abs(diff) > 0.1 & PBSvsWT > 0 & Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT >= 0 & PBSvsWT > Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT ~ "rescue",
                             abs(diff) > 0.1 & PBSvsWT < 0 & Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT <= 0 & PBSvsWT < Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT ~ "rescue",
                             abs(diff) > 0.1 & PBSvsWT > 0 & Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT >= 0 & PBSvsWT < Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT ~ "reversal",
                             abs(diff) > 0.1 & PBSvsWT < 0 & Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT <= 0 & PBSvsWT > Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT ~ "reversal",
                             abs(diff) > 0.1 & PBSvsWT > 0 & 
                               PBSvsWT > (Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT) & 
                               (Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT >= -0.2 & Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT <= 0.2) ~ "rescue",
                             abs(diff) > 0.1 & PBSvsWT > 0 &  
                               PBSvsWT > (Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT) & 
                               (Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT < -0.2 | Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT > 0.2) ~ "reversal",
                             abs(diff) > 0.1 & PBSvsWT < 0 & 
                               (PBSvsWT) < Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT & 
                               (Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT >= -0.2 & Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT <= 0.2) ~ "rescue",
                             abs(diff) > 0.1 & PBSvsWT < 0 &  
                               (PBSvsWT) < Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT & 
                               (Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT < -0.2 | Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT > 0.2) ~ "reversal"))
# Replace NA 
df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_sig_wd_delta <- 
  df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_sig_wd_delta %>% 
  replace_na(list(SJ_type = "diff_inf_0.1"))
# long format to plot
df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_SJmore1_sig_lg_delta <- pivot_longer(df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_sig_wd_delta,
                                                                                                              cols = c(PBSvsWT, Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT),
                                                                                                              names_to = "comparison",
                                                                                                              values_to = "IncLevelDifference")
df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_SJmore1_sig_lg_delta <- 
  df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_SJmore1_sig_lg_delta %>% 
  mutate(comparison = factor(comparison, levels = c("PBSvsWT", "Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT")))
# Add Ttn column
df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_SJmore1_sig_lg_delta_ttn <- 
  df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_SJmore1_sig_lg_delta %>% 
  mutate("Ttn" = geneSymbol == "Ttn")
# Add color variable
df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_SJmore1_sig_lg_delta_ttn <- 
  df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_SJmore1_sig_lg_delta_ttn %>% 
  mutate(plot_color = if_else(Ttn == TRUE, true = "Ttn", false = as.character(SJ_type)))
# Factor order
df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_SJmore1_sig_lg_delta_ttn <- 
  df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_SJmore1_sig_lg_delta_ttn %>% 
  mutate(SJ_type = factor(SJ_type, levels = c("rescue", "reversal", "diff_inf_0.1")))
# Plot
pdf("rmats_directory/VF_P635L_after_base_editing_dotlineplot_allSJs_common_sigin_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_deltadeltaPSI_facet.pdf", 
    width = 4, height = 6)
ggplot(df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_SJmore1_sig_lg_delta_ttn,
       aes(x = comparison, y = IncLevelDifference, group = SJ_ID, color = plot_color)) +
  geom_point(alpha = 0.6, aes(size = plot_color)) +
  geom_line(alpha = 0.6, aes(size = plot_color)) + 
  facet_wrap(vars(SJ_type), strip.position="right", ncol = 1, labeller = as_labeller(c("rescue" = "Rescue",
                                                                                       "reversal" = "Mis-spliced", 
                                                                                       "diff_inf_0.1" = "Unchanged"))) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  scale_color_manual(name = NULL, 
                     values = c("diff_inf_0.1" = "gray58", "rescue" = "red", "reversal" = "gray30", "Ttn" = "blue"),
                     breaks = c("Ttn")) +
  scale_size_manual(guide = 'none', values = c(0.2, 0.2, 0.2, 0.8)) +
  ggtitle("P635L") +
  ylab("Delta delta PSI") +
  scale_x_discrete(name = element_blank(), 
                   labels = c("PBS", "SpRY")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  annotate("rect", xmin = 0, xmax = 3, ymin = -0.2, ymax = 0.2, 
           alpha = .2) +
  guides(colour = guide_legend(override.aes = list(alpha = 1,
                                                   linetype = NULL,
                                                   size = 2))) 
dev.off()



#### Number of events in each category ####
df_nbevents_per_class <- bind_rows(data.frame("mutant" = "P635L", 
                     "comparison" = "PBS - SpRy",
                     table(unique(df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_SJmore1_sig_lg_delta_ttn[,c("SJ_ID", "SJ_type")])$SJ_type)),
          data.frame("mutant" = "P635L", 
                     "comparison" = "PBS - Abe8e",
                     table(unique(df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_SJmore1_sig_lg_delta_ttn[,c("SJ_ID", "SJ_type")])$SJ_type)), 
          data.frame("mutant" = "R636Q", 
                     "comparison" = "PBS - Abe8e",
                     table(unique(df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_sig_lg_delta_ttn[,c("SJ_ID", "SJ_type")])$SJ_type)))
# Modify names
df_nbevents_per_class <- mutate(df_nbevents_per_class, 
       Var1 = case_when(Var1 == "reversal" ~ "Mis-spliced",
                        Var1 == "rescue" ~ "Rescue",
                        Var1 == "diff_inf_0.1" ~ "Unchanged"))
# wide format
df_nbevents_per_class_wd <- 
  pivot_wider(df_nbevents_per_class, 
              id_cols = c("mutant", "comparison"), 
              names_from = "Var1", 
              values_from = "Freq")
write.table(df_nbevents_per_class_wd,
            file = "rmats_directory/table_nbevents_per_class.txt",
            col.names = TRUE,
            row.names = FALSE)

#### Extract events that correspond to RBM20 targets ####
# Read list from Markus
v_rbm20_genes <- read_xlsx("results/RBM20_splice_targets_Hoogenhof.xlsx", col_names = FALSE)
# Format gene names
library(stringr)
v_rbm20_genes <- str_to_title(as.vector(v_rbm20_genes$`...1`))
# Identify genes in the rmats analysis
df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_rbm20target <- 
  df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1 %>% 
  dplyr::filter(geneSymbol %in% v_rbm20_genes)
df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_rbm20target <- 
  df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1 %>% 
  dplyr::filter(geneSymbol %in% v_rbm20_genes)
df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_rbm20target <-
  df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1 %>% 
  dplyr::filter(geneSymbol %in% v_rbm20_genes)
# Merge
df_merge_rbm20target <-
  dplyr::bind_rows(list(R636Q_Abe8e=df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_rbm20target, 
                      P635L_Abe8e=df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_rbm20target,
                      P635L_SpRY=df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_rbm20target), 
                 .id = 'source')
# Identify events that overlap
df_merge_rbm20target_SJID <- df_merge_rbm20target %>% 
  group_by(geneSymbol, SJ_name, longExonStart_0base, longExonEnd,
           shortES, shortEE, flankingES, flankingEE,
           X1stExonStart_0base, X1stExonEnd, X2ndExonStart_0base, X2ndExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE,
           riExonEnd, exonStart_0base, exonEnd) %>% 
  mutate(merge_SJ_ID = cur_group_id())
# Filter for SJ events found more than once
df_merge_rbm20target_SJID_in6 <- df_merge_rbm20target_SJID %>% 
  group_by(merge_SJ_ID) %>% 
  dplyr::filter(n() == 6)
# Significant
df_merge_rbm20target_SJID_in6_sig <- df_merge_rbm20target_SJID_in6 %>% 
  group_by(merge_SJ_ID) %>% 
  dplyr::filter((all(FDR[comparison == "PBSvsWT"] < 0.01) | all(FDR[comparison == "PBSvsWT"] == 0)) &
                  (all(IncLevelDifference[comparison == "PBSvsWT"] < -0.1) | 
                     all(IncLevelDifference[comparison == "PBSvsWT"] > 0.1)))
write.table(df_merge_rbm20target_SJID_in6_sig, 
            file = "rmats_directory/df_merge_rbm20target_SJID_in6_sig.txt",
            row.names = FALSE)
# Add a column to do the wide format
df_merge_rbm20target_SJID_in6_sig <- df_merge_rbm20target_SJID_in6_sig %>% 
  mutate(comparison_simp = ifelse(comparison == "PBSvsWT", "PBSvsWT", "editedvsWT"))
# Wide format for Markus
df_merge_rbm20target_SJID_in6_sig_wd <- df_merge_rbm20target_SJID_in6_sig %>% 
  dplyr::select(source, merge_SJ_ID, geneSymbol, SJ_ID, IncLevelDifference, comparison_simp) %>% 
  pivot_wider(names_from = comparison_simp, values_from = IncLevelDifference)
# SJ events classification
# Difference
df_merge_rbm20target_SJID_in6_sig_wd_delta <- 
  df_merge_rbm20target_SJID_in6_sig_wd %>% 
  mutate(diff = PBSvsWT - editedvsWT) 
# Classification
df_merge_rbm20target_SJID_in6_sig_wd_delta <- 
  df_merge_rbm20target_SJID_in6_sig_wd_delta %>% 
  mutate(SJ_type = case_when(abs(diff) > 0.1 & PBSvsWT > 0 & editedvsWT >= 0 & PBSvsWT > editedvsWT ~ "rescue",
                             abs(diff) > 0.1 & PBSvsWT < 0 & editedvsWT <= 0 & PBSvsWT < editedvsWT ~ "rescue",
                             abs(diff) > 0.1 & PBSvsWT > 0 & editedvsWT >= 0 & PBSvsWT < editedvsWT ~ "reversal",
                             abs(diff) > 0.1 & PBSvsWT < 0 & editedvsWT <= 0 & PBSvsWT > editedvsWT ~ "reversal",
                             abs(diff) > 0.1 & PBSvsWT > 0 & 
                               PBSvsWT > (editedvsWT) & 
                               (editedvsWT >= -0.2 & editedvsWT <= 0.2) ~ "rescue",
                             abs(diff) > 0.1 & PBSvsWT > 0 &  
                               PBSvsWT > (editedvsWT) & 
                               (editedvsWT < -0.2 | editedvsWT > 0.2) ~ "reversal",
                             abs(diff) > 0.1 & PBSvsWT < 0 & 
                               (PBSvsWT) < editedvsWT & 
                               (editedvsWT >= -0.2 & editedvsWT <= 0.2) ~ "rescue",
                             abs(diff) > 0.1 & PBSvsWT < 0 &  
                               (PBSvsWT) < editedvsWT & 
                               (editedvsWT < -0.2 | editedvsWT > 0.2) ~ "reversal"))
# Replace NA 
df_merge_rbm20target_SJID_in6_sig_wd_delta <- 
  df_merge_rbm20target_SJID_in6_sig_wd_delta %>% 
  replace_na(list(SJ_type = "diff_inf_0.1"))
# long format to plot
df_merge_rbm20target_SJID_in6_sig_wd_delta_lg <- pivot_longer(df_merge_rbm20target_SJID_in6_sig_wd_delta,
                                                              cols = c(PBSvsWT, editedvsWT),
                                                              names_to = "comparison",
                                                              values_to = "IncLevelDifference")

### Markus approach ###
# Extract all genes from the RBM20 targets from PBS vs WT comparisons
v_R636Q_after_base_editing_PBS_SJmore1_rbm20target_genesymbols <- 
  df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1 %>% 
  dplyr::filter(comparison == "PBSvsWT") %>% 
  dplyr::filter(geneSymbol %in% v_rbm20_genes) %>%
  pull(geneSymbol) %>%
  unique()
v_P635L_after_base_editing_PBS_Abe8e_SJmore1_rbm20target_genesymbols <- 
  df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1 %>% 
  dplyr::filter(comparison == "PBSvsWT") %>% 
  dplyr::filter(geneSymbol %in% v_rbm20_genes) %>%
  pull(geneSymbol) %>%
  unique()
v_P635L_after_base_editing_PBS_SpRY_SJmore1_rbm20target_genesymbols <- 
  df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1 %>% 
  dplyr::filter(comparison == "PBSvsWT") %>% 
  dplyr::filter(geneSymbol %in% v_rbm20_genes) %>%
  pull(geneSymbol) %>%
  unique()
# Filter corresponding events and events that overlap
df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_rbm20target <- 
  df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1 %>% 
  dplyr::filter(geneSymbol %in% v_R636Q_after_base_editing_PBS_SJmore1_rbm20target_genesymbols) %>% 
  group_by(geneSymbol, SJ_name, longExonStart_0base, longExonEnd,
           shortES, shortEE, flankingES, flankingEE,
           X1stExonStart_0base, X1stExonEnd, X2ndExonStart_0base, X2ndExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE,
           riExonEnd, exonStart_0base, exonEnd) %>% 
  mutate(merge_SJ_ID = cur_group_id())
df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_rbm20target <- 
  df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1 %>% 
  dplyr::filter(geneSymbol %in% v_P635L_after_base_editing_PBS_Abe8e_SJmore1_rbm20target_genesymbols) %>% 
  group_by(geneSymbol, SJ_name, longExonStart_0base, longExonEnd,
           shortES, shortEE, flankingES, flankingEE,
           X1stExonStart_0base, X1stExonEnd, X2ndExonStart_0base, X2ndExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE,
           riExonEnd, exonStart_0base, exonEnd) %>% 
  mutate(merge_SJ_ID = cur_group_id())
mutate(merge_SJ_ID = cur_group_id())
df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_rbm20target <- 
  df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1 %>% 
  dplyr::filter(geneSymbol %in% v_P635L_after_base_editing_PBS_SpRY_SJmore1_rbm20target_genesymbols) %>% 
  group_by(geneSymbol, SJ_name, longExonStart_0base, longExonEnd,
           shortES, shortEE, flankingES, flankingEE,
           X1stExonStart_0base, X1stExonEnd, X2ndExonStart_0base, X2ndExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE,
           riExonEnd, exonStart_0base, exonEnd) %>% 
  mutate(merge_SJ_ID = cur_group_id())
# Compute mean per event
df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_rbm20target_wmean <- 
  df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_rbm20target %>% 
  rowwise() %>% 
  mutate(avg_IncLevel1 = mean(as.numeric(unlist(strsplit(IncLevel1, ","))), na.rm = TRUE),
         avg_IncLevel2 = mean(as.numeric(unlist(strsplit(IncLevel2, ","))), na.rm = TRUE))
df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_rbm20target_wmean <- 
  df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_rbm20target %>% 
  rowwise() %>% 
  mutate(avg_IncLevel1 = mean(as.numeric(unlist(strsplit(IncLevel1, ","))), na.rm = TRUE),
         avg_IncLevel2 = mean(as.numeric(unlist(strsplit(IncLevel2, ","))), na.rm = TRUE))
df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_rbm20target_wmean <- 
  df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_rbm20target %>% 
  rowwise() %>% 
  mutate(avg_IncLevel1 = mean(as.numeric(unlist(strsplit(IncLevel1, ","))), na.rm = TRUE),
         avg_IncLevel2 = mean(as.numeric(unlist(strsplit(IncLevel2, ","))), na.rm = TRUE))
# Plot
### P635L SpRy ###
# Prepare data to plot
# test <- as.data.frame(df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_rbm20target_wmean) %>% 
#   dplyr::select(geneSymbol, comparison, FDR, mutant, merge_SJ_ID, avg_IncLevel1, avg_IncLevel2) %>%
#   pivot_longer(cols = c("avg_IncLevel1", "avg_IncLevel2"), names_to = "PSI", values_to = "mean")
df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_rbm20target_wmean_w <- as.data.frame(df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_rbm20target_wmean) %>% 
  dplyr::select(geneSymbol, comparison, FDR, mutant, merge_SJ_ID, avg_IncLevel1, avg_IncLevel2) %>%
  pivot_wider(names_from = comparison,
              names_glue = "{comparison}_{.value}",
              values_from = c(avg_IncLevel1, avg_IncLevel2, FDR))
# All events
pdf("rmats_directory/heatmap_rmats_P635L_rbm20target_Hoogenhof_SpRy.pdf", width = 3, height = 15)
pheatmap(df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_rbm20target_wmean_w[,c("PBSvsWT_avg_IncLevel2",
                 "Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_avg_IncLevel1",
                 "PBSvsWT_avg_IncLevel1")],
         color=colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         labels_row = df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_rbm20target_wmean_w$geneSymbol,
         labels_col = c("WT", "SpRy", "PBS"),
         na_col = "gray58",
         fontsize = 5, 
         main = "P635L - Hoogenhof list")
dev.off()
# Events significant in PBS vs WT
df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_rbm20target_wmean_w_sig <- df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_rbm20target_wmean_w %>% 
  dplyr::filter(PBSvsWT_FDR < 0.01 | PBSvsWT_FDR == 0) 
pdf("rmats_directory/heatmap_rmats_P635L_rbm20target_Hoogenhof_SpRy_siginPBSvsWT.pdf", width = 3, height = 6)
pheatmap(df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_rbm20target_wmean_w_sig[,c("PBSvsWT_avg_IncLevel2",
                 "Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_avg_IncLevel1",
                 "PBSvsWT_avg_IncLevel1")],
         color=colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         labels_row = df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_rbm20target_wmean_w$geneSymbol,
         labels_col = c("WT", "SpRy", "PBS"),
         na_col = "gray58",
         fontsize = 5, 
         main = "P635L - Hoogenhof list - significant events in PBS vs WT")
dev.off()
# Events significant in both PBS vs WT and SpRY vs WT
df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_rbm20target_wmean_w_sig_both <- df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_rbm20target_wmean_w %>% 
  dplyr::filter((PBSvsWT_FDR < 0.01 | PBSvsWT_FDR == 0) & 
                  (Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_FDR < 0.01 | Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_FDR == 0)) 
pdf("rmats_directory/heatmap_rmats_P635L_rbm20target_Hoogenhof_SpRy_siginPBSvsWTandSpRYvsWT.pdf", 
    width = 4, height = 3)
pheatmap(df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_rbm20target_wmean_w_sig_both[,c("PBSvsWT_avg_IncLevel2",
                     "Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_avg_IncLevel1",
                     "PBSvsWT_avg_IncLevel1")],
         color=colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         labels_row = df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_rbm20target_wmean_w$geneSymbol,
         labels_col = c("WT", "SpRy", "PBS"),
         na_col = "gray58",
         fontsize = 5, 
         main = "P635L - Hoogenhof list - significant events in PBS vs WT and SpRY vs WT")
dev.off()

### P635L Abe8e ###
# Prepare data to plot
df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_rbm20target_wmean_w <- as.data.frame(df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_rbm20target_wmean) %>% 
  dplyr::select(geneSymbol, comparison, FDR, mutant, merge_SJ_ID, avg_IncLevel1, avg_IncLevel2) %>%
  pivot_wider(names_from = comparison,
              names_glue = "{comparison}_{.value}",
              values_from = c(avg_IncLevel1, avg_IncLevel2, FDR))
# All events
pdf("rmats_directory/heatmap_rmats_P635L_rbm20target_Hoogenhof_Abe8e.pdf", width = 3, height = 15)
pheatmap(df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_rbm20target_wmean_w[,c("PBSvsWT_avg_IncLevel2",
                                                                                                            "Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_avg_IncLevel1",
                                                                                                            "PBSvsWT_avg_IncLevel1")],
         color=colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         labels_row = df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_rbm20target_wmean_w$geneSymbol,
         labels_col = c("WT", "Abe8e", "PBS"),
         na_col = "gray58",
         fontsize = 5, 
         main = "P635L - Hoogenhof list")
dev.off()
# Events significant in PBS vs WT
df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_rbm20target_wmean_w_sig <- df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_rbm20target_wmean_w %>% 
  dplyr::filter(PBSvsWT_FDR < 0.01 | PBSvsWT_FDR == 0) 
pdf("rmats_directory/heatmap_rmats_P635L_rbm20target_Hoogenhof_Abe8e_siginPBSvsWT.pdf", width = 3, height = 6)
pheatmap(df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_rbm20target_wmean_w_sig[,c("PBSvsWT_avg_IncLevel2",
                                                                                                                "Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_avg_IncLevel1",
                                                                                                                "PBSvsWT_avg_IncLevel1")],
         color=colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         labels_row = df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_rbm20target_wmean_w$geneSymbol,
         labels_col = c("WT", "Abe8e", "PBS"),
         na_col = "gray58",
         fontsize = 5, 
         main = "P635L - Hoogenhof list - significant events in PBS vs WT")
dev.off()
# Events significant in both PBS vs WT and Abe8e vs WT
df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_rbm20target_wmean_w_sig_both <- df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_rbm20target_wmean_w %>% 
  dplyr::filter((PBSvsWT_FDR < 0.01 | PBSvsWT_FDR == 0) & 
                  (Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_FDR < 0.01 | Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_FDR == 0)) 
pdf("rmats_directory/heatmap_rmats_P635L_rbm20target_Hoogenhof_Abe8e_siginPBSvsWTandAbe8evsWT.pdf", 
    width = 4, height = 3)
pheatmap(df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_rbm20target_wmean_w_sig_both[,c("PBSvsWT_avg_IncLevel2",
                                                                                                                     "Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_avg_IncLevel1",
                                                                                                                     "PBSvsWT_avg_IncLevel1")],
         color=colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         labels_row = df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_rbm20target_wmean_w$geneSymbol,
         labels_col = c("WT", "Abe8e", "PBS"),
         na_col = "gray58",
         fontsize = 5, 
         main = "P635L - Hoogenhof list - significant events in PBS vs WT and Abe8e vs WT")
dev.off()

### R636Q Abe8e ###
# Prepare data to plot
df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_rbm20target_wmean_w <- as.data.frame(df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_rbm20target_wmean) %>% 
  dplyr::select(geneSymbol, comparison, FDR, mutant, merge_SJ_ID, avg_IncLevel1, avg_IncLevel2) %>%
  pivot_wider(names_from = comparison,
              names_glue = "{comparison}_{.value}",
              values_from = c(avg_IncLevel1, avg_IncLevel2, FDR))
# All events
pdf("rmats_directory/heatmap_rmats_R636Q_rbm20target_Hoogenhof_Abe8e.pdf", width = 3, height = 15)
pheatmap(df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_rbm20target_wmean_w[,c("PBSvsWT_avg_IncLevel2",
                                                                                                             "Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_avg_IncLevel1",
                                                                                                             "PBSvsWT_avg_IncLevel1")],
         color=colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         labels_row = df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_rbm20target_wmean_w$geneSymbol,
         labels_col = c("WT", "Abe8e", "PBS"),
         na_col = "gray58",
         fontsize = 5, 
         main = "R636Q - Hoogenhof list")
dev.off()
# Events significant in PBS vs WT
df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_rbm20target_wmean_w_sig <- df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_rbm20target_wmean_w %>% 
  dplyr::filter(PBSvsWT_FDR < 0.01 | PBSvsWT_FDR == 0) 
pdf("rmats_directory/heatmap_rmats_R636Q_rbm20target_Hoogenhof_Abe8e_siginPBSvsWT.pdf", width = 3, height = 6)
pheatmap(df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_rbm20target_wmean_w_sig[,c("PBSvsWT_avg_IncLevel2",
                                                                                                                 "Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_avg_IncLevel1",
                                                                                                                 "PBSvsWT_avg_IncLevel1")],
         color=colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         labels_row = df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_rbm20target_wmean_w$geneSymbol,
         labels_col = c("WT", "Abe8e", "PBS"),
         na_col = "gray58",
         fontsize = 5, 
         main = "R636Q - Hoogenhof list - significant events in PBS vs WT")
dev.off()
# Events significant in both PBS vs WT and Abe8e vs WT
df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_rbm20target_wmean_w_sig_both <- df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_rbm20target_wmean_w %>% 
  dplyr::filter((PBSvsWT_FDR < 0.01 | PBSvsWT_FDR == 0) & 
                  (Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_FDR < 0.01 | Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_FDR == 0)) 
pdf("rmats_directory/heatmap_rmats_R636Q_rbm20target_Hoogenhof_Abe8e_siginPBSvsWTandAbe8evsWT.pdf", 
    width = 4, height = 3)
pheatmap(df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_rbm20target_wmean_w_sig_both[,c("PBSvsWT_avg_IncLevel2",
                                                                                                                      "Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_avg_IncLevel1",
                                                                                                                      "PBSvsWT_avg_IncLevel1")],
         color=colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         labels_row = df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_rbm20target_wmean_w$geneSymbol,
         labels_col = c("WT", "Abe8e", "PBS"),
         na_col = "gray58",
         fontsize = 5, 
         main = "R636Q - Hoogenhof list - significant events in PBS vs WT and Abe8e vs WT")
dev.off()

# DEE overlap per Markus from the DESeq2 analysis
v_DEG_DESeq <- as.vector(read.table("rmats_directory/Final_List_DEEs_overlap_P635L_HOM_R636Q_HOM.txt", 
           header = FALSE))$V1
# Extract all genes in PBS vs WT comparisons
v_R636Q_after_base_editing_PBS_SJmore1_DEG_DESeq_genesymbols <- 
  df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1 %>% 
  dplyr::filter(comparison == "PBSvsWT") %>% 
  dplyr::filter(geneSymbol %in% v_DEG_DESeq) %>%
  pull(geneSymbol) %>%
  unique()
v_P635L_after_base_editing_PBS_Abe8e_SJmore1_DEG_DESeq_genesymbols <- 
  df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1 %>% 
  dplyr::filter(comparison == "PBSvsWT") %>% 
  dplyr::filter(geneSymbol %in% v_DEG_DESeq) %>%
  pull(geneSymbol) %>%
  unique()
v_P635L_after_base_editing_PBS_SpRY_SJmore1_DEG_DESeq_genesymbols <- 
  df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1 %>% 
  dplyr::filter(comparison == "PBSvsWT") %>% 
  dplyr::filter(geneSymbol %in% v_DEG_DESeq) %>%
  pull(geneSymbol) %>%
  unique()
# Filter corresponding events and events that overlap
df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_DEG_DESeq <- 
  df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1 %>% 
  dplyr::filter(geneSymbol %in% v_R636Q_after_base_editing_PBS_SJmore1_DEG_DESeq_genesymbols) %>% 
  group_by(geneSymbol, SJ_name, longExonStart_0base, longExonEnd,
           shortES, shortEE, flankingES, flankingEE,
           X1stExonStart_0base, X1stExonEnd, X2ndExonStart_0base, X2ndExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE,
           riExonEnd, exonStart_0base, exonEnd) %>% 
  mutate(merge_SJ_ID = cur_group_id())
df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_DEG_DESeq <- 
  df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1 %>% 
  dplyr::filter(geneSymbol %in% v_P635L_after_base_editing_PBS_Abe8e_SJmore1_DEG_DESeq_genesymbols) %>% 
  group_by(geneSymbol, SJ_name, longExonStart_0base, longExonEnd,
           shortES, shortEE, flankingES, flankingEE,
           X1stExonStart_0base, X1stExonEnd, X2ndExonStart_0base, X2ndExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE,
           riExonEnd, exonStart_0base, exonEnd) %>% 
  mutate(merge_SJ_ID = cur_group_id())
mutate(merge_SJ_ID = cur_group_id())
df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_DEG_DESeq <- 
  df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1 %>% 
  dplyr::filter(geneSymbol %in% v_P635L_after_base_editing_PBS_SpRY_SJmore1_DEG_DESeq_genesymbols) %>% 
  group_by(geneSymbol, SJ_name, longExonStart_0base, longExonEnd,
           shortES, shortEE, flankingES, flankingEE,
           X1stExonStart_0base, X1stExonEnd, X2ndExonStart_0base, X2ndExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE,
           riExonEnd, exonStart_0base, exonEnd) %>% 
  mutate(merge_SJ_ID = cur_group_id())
# Compute mean per event
df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_DEG_DESeq_wmean <- 
  df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_DEG_DESeq %>% 
  rowwise() %>% 
  mutate(avg_IncLevel1 = mean(as.numeric(unlist(strsplit(IncLevel1, ","))), na.rm = TRUE),
         avg_IncLevel2 = mean(as.numeric(unlist(strsplit(IncLevel2, ","))), na.rm = TRUE))
df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_DEG_DESeq_wmean <- 
  df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_DEG_DESeq %>% 
  rowwise() %>% 
  mutate(avg_IncLevel1 = mean(as.numeric(unlist(strsplit(IncLevel1, ","))), na.rm = TRUE),
         avg_IncLevel2 = mean(as.numeric(unlist(strsplit(IncLevel2, ","))), na.rm = TRUE))
df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_DEG_DESeq_wmean <- 
  df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_DEG_DESeq %>% 
  rowwise() %>% 
  mutate(avg_IncLevel1 = mean(as.numeric(unlist(strsplit(IncLevel1, ","))), na.rm = TRUE),
         avg_IncLevel2 = mean(as.numeric(unlist(strsplit(IncLevel2, ","))), na.rm = TRUE))
# Plot
### P635L SpRy ###
# Prepare data to plot
df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_DEG_DESeq_wmean_w <- as.data.frame(df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_DEG_DESeq_wmean) %>% 
  dplyr::select(geneSymbol, comparison, FDR, mutant, merge_SJ_ID, avg_IncLevel1, avg_IncLevel2) %>%
  pivot_wider(names_from = comparison,
              names_glue = "{comparison}_{.value}",
              values_from = c(avg_IncLevel1, avg_IncLevel2, FDR))
# All events
pdf("rmats_directory/heatmap_rmats_P635L_DEG_DESeq_SpRy.pdf", width = 3, height = 15)
pheatmap(df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_DEG_DESeq_wmean_w[,c("PBSvsWT_avg_IncLevel2",
                                                                                                            "Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_avg_IncLevel1",
                                                                                                            "PBSvsWT_avg_IncLevel1")],
         color=colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         labels_row = df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_DEG_DESeq_wmean_w$geneSymbol,
         labels_col = c("WT", "SpRy", "PBS"),
         na_col = "gray58",
         fontsize = 5, 
         main = "P635L - DEG DESeq list")
dev.off()
# Events significant in PBS vs WT
df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_DEG_DESeq_wmean_w_sig <- df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_DEG_DESeq_wmean_w %>% 
  dplyr::filter(PBSvsWT_FDR < 0.01 | PBSvsWT_FDR == 0) 
pdf("rmats_directory/heatmap_rmats_P635L_DEG_DESeq_SpRy_siginPBSvsWT.pdf", width = 3, height = 6)
pheatmap(df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_DEG_DESeq_wmean_w_sig[,c("PBSvsWT_avg_IncLevel2",
                                                                                                                "Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_avg_IncLevel1",
                                                                                                                "PBSvsWT_avg_IncLevel1")],
         color=colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         labels_row = df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_DEG_DESeq_wmean_w$geneSymbol,
         labels_col = c("WT", "SpRy", "PBS"),
         na_col = "gray58",
         fontsize = 5, 
         main = "P635L - DEG DESeq list - significant events in PBS vs WT")
dev.off()
# Events significant in both PBS vs WT and SpRY vs WT
df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_DEG_DESeq_wmean_w_sig_both <- df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_DEG_DESeq_wmean_w %>% 
  dplyr::filter((PBSvsWT_FDR < 0.01 | PBSvsWT_FDR == 0) & 
                  (Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_FDR < 0.01 | Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_FDR == 0)) 
pdf("rmats_directory/heatmap_rmats_P635L_DEG_DESeq_SpRy_siginPBSvsWTandSpRYvsWT.pdf", 
    width = 4, height = 3)
pheatmap(df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_DEG_DESeq_wmean_w_sig_both[,c("PBSvsWT_avg_IncLevel2",
                                                                                                                     "Nterm_SpRY_and_Cterm_SpRY_gRNA5vsWT_avg_IncLevel1",
                                                                                                                     "PBSvsWT_avg_IncLevel1")],
         color=colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         labels_row = df_P635L_after_base_editing_PBS_and_Nterm_SpRY_and_Cterm_SpRY_gRNA5_SJmore1_DEG_DESeq_wmean_w$geneSymbol,
         labels_col = c("WT", "SpRy", "PBS"),
         na_col = "gray58",
         fontsize = 5, 
         main = "P635L - DEG DESeq list - significant events in PBS vs WT and SpRY vs WT")
dev.off()

### P635L Abe8e ###
# Prepare data to plot
df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_DEG_DESeq_wmean_w <- as.data.frame(df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_DEG_DESeq_wmean) %>% 
  dplyr::select(geneSymbol, comparison, FDR, mutant, merge_SJ_ID, avg_IncLevel1, avg_IncLevel2) %>%
  pivot_wider(names_from = comparison,
              names_glue = "{comparison}_{.value}",
              values_from = c(avg_IncLevel1, avg_IncLevel2, FDR))
# All events
pdf("rmats_directory/heatmap_rmats_P635L_DEG_DESeq_Abe8e.pdf", width = 3, height = 15)
pheatmap(df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_DEG_DESeq_wmean_w[,c("PBSvsWT_avg_IncLevel2",
                                                                                                             "Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_avg_IncLevel1",
                                                                                                             "PBSvsWT_avg_IncLevel1")],
         color=colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         labels_row = df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_DEG_DESeq_wmean_w$geneSymbol,
         labels_col = c("WT", "Abe8e", "PBS"),
         na_col = "gray58",
         fontsize = 5, 
         main = "P635L - DEG DESeq list")
dev.off()
# Events significant in PBS vs WT
df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_DEG_DESeq_wmean_w_sig <- df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_DEG_DESeq_wmean_w %>% 
  dplyr::filter(PBSvsWT_FDR < 0.01 | PBSvsWT_FDR == 0) 
pdf("rmats_directory/heatmap_rmats_P635L_DEG_DESeq_Abe8e_siginPBSvsWT.pdf", width = 3, height = 6)
pheatmap(df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_DEG_DESeq_wmean_w_sig[,c("PBSvsWT_avg_IncLevel2",
                                                                                                                 "Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_avg_IncLevel1",
                                                                                                                 "PBSvsWT_avg_IncLevel1")],
         color=colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         labels_row = df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_DEG_DESeq_wmean_w$geneSymbol,
         labels_col = c("WT", "Abe8e", "PBS"),
         na_col = "gray58",
         fontsize = 5, 
         main = "P635L - DEG DESeq list - significant events in PBS vs WT")
dev.off()
# Events significant in both PBS vs WT and Abe8e vs WT
df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_DEG_DESeq_wmean_w_sig_both <- df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_DEG_DESeq_wmean_w %>% 
  dplyr::filter((PBSvsWT_FDR < 0.01 | PBSvsWT_FDR == 0) & 
                  (Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_FDR < 0.01 | Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_FDR == 0)) 
pdf("rmats_directory/heatmap_rmats_P635L_DEG_DESeq_Abe8e_siginPBSvsWTandAbe8evsWT.pdf", 
    width = 4, height = 3)
pheatmap(df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_DEG_DESeq_wmean_w_sig_both[,c("PBSvsWT_avg_IncLevel2",
                                                                                                                      "Nterm_NRTH_Abe8e_and_Cterm_gRNA5vsWT_avg_IncLevel1",
                                                                                                                      "PBSvsWT_avg_IncLevel1")],
         color=colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         labels_row = df_P635L_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_gRNA5_SJmore1_DEG_DESeq_wmean_w$geneSymbol,
         labels_col = c("WT", "Abe8e", "PBS"),
         na_col = "gray58",
         fontsize = 5, 
         main = "P635L - DEG DESeq list - significant events in PBS vs WT and Abe8e vs WT")
dev.off()

### R636Q Abe8e ###
# Prepare data to plot
df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_DEG_DESeq_wmean_w <- as.data.frame(df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_DEG_DESeq_wmean) %>% 
  dplyr::select(geneSymbol, comparison, FDR, mutant, merge_SJ_ID, avg_IncLevel1, avg_IncLevel2) %>%
  pivot_wider(names_from = comparison,
              names_glue = "{comparison}_{.value}",
              values_from = c(avg_IncLevel1, avg_IncLevel2, FDR))
# All events
pdf("rmats_directory/heatmap_rmats_R636Q_DEG_DESeq_Abe8e.pdf", width = 3, height = 15)
pheatmap(df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_DEG_DESeq_wmean_w[,c("PBSvsWT_avg_IncLevel2",
                                                                                                                  "Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_avg_IncLevel1",
                                                                                                                  "PBSvsWT_avg_IncLevel1")],
         color=colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         labels_row = df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_DEG_DESeq_wmean_w$geneSymbol,
         labels_col = c("WT", "Abe8e", "PBS"),
         na_col = "gray58",
         fontsize = 5, 
         main = "R636Q - DEG DESeq list")
dev.off()
# Events significant in PBS vs WT
df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_DEG_DESeq_wmean_w_sig <- df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_DEG_DESeq_wmean_w %>% 
  dplyr::filter(PBSvsWT_FDR < 0.01 | PBSvsWT_FDR == 0) 
pdf("rmats_directory/heatmap_rmats_R636Q_DEG_DESeq_Abe8e_siginPBSvsWT.pdf", width = 3, height = 6)
pheatmap(df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_DEG_DESeq_wmean_w_sig[,c("PBSvsWT_avg_IncLevel2",
                                                                                                                      "Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_avg_IncLevel1",
                                                                                                                      "PBSvsWT_avg_IncLevel1")],
         color=colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         labels_row = df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_DEG_DESeq_wmean_w$geneSymbol,
         labels_col = c("WT", "Abe8e", "PBS"),
         na_col = "gray58",
         fontsize = 5, 
         main = "R636Q - DEG DESeq list - significant events in PBS vs WT")
dev.off()
# Events significant in both PBS vs WT and Abe8e vs WT
df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_DEG_DESeq_wmean_w_sig_both <- df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_DEG_DESeq_wmean_w %>% 
  dplyr::filter((PBSvsWT_FDR < 0.01 | PBSvsWT_FDR == 0) & 
                  (Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_FDR < 0.01 | Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_FDR == 0)) 
pdf("rmats_directory/heatmap_rmats_R636Q_DEG_DESeq_Abe8e_siginPBSvsWTandAbe8evsWT.pdf", 
    width = 4, height = 3)
pheatmap(df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_DEG_DESeq_wmean_w_sig_both[,c("PBSvsWT_avg_IncLevel2",
                                                                                                                           "Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_avg_IncLevel1",
                                                                                                                           "PBSvsWT_avg_IncLevel1")],
         color=colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         labels_row = df_R636Q_after_base_editing_PBS_and_Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6_SJmore1_DEG_DESeq_wmean_w$geneSymbol,
         labels_col = c("WT", "Abe8e", "PBS"),
         na_col = "gray58",
         fontsize = 5, 
         main = "R636Q - DEG DESeq list - significant events in PBS vs WT and Abe8e vs WT")
dev.off()
