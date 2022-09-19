# Create heatmap of gene expression

# Packages
library(pheatmap)
library(tidyverse)
library(ggplot2)

# Load the data
# After base editing
featurecounts_P635L_after_base_editing <- 
  readRDS("results/robject_featurecounts_P635L_after_base_editing.rds")
# No base editing
# featurecounts_P635L <- 
#   readRDS("results/robject_featurecounts_P635L.rds")

# Add metadata
# Sample table
# sample_table_noedit <- read.table("/g/steinmetz/project/dcm_lgreads/mouse_bulkRNA/config/samples_metadata.txt",
#                                   header = TRUE)
sample_table_wedit <- read.table("/g/steinmetz/project/dcm_lgreads/mouse_bulkRNA/config/samples_afterbaseediting_metadata.txt",
                                 header = TRUE)
# Add info
# Simplify sample table
# sample_table_noedit$sample <- paste("results/STAR", 
#                                     sample_table_noedit$mutant, 
#                                     sample_table_noedit$run, 
#                                     "Aligned.sortedByCoord.out.bam", sep = "/")
sample_table_wedit$sample <- paste("results/STAR", 
                            sample_table_wedit$mutant, 
                            sample_table_wedit$run, 
                            "Aligned.sortedByCoord.out.bam", sep = "/")
# Filter
# sample_table_noedit_P635L <- dplyr::filter(sample_table_noedit, mutant == "P635L")
sample_table_wedit_P635L <- dplyr::filter(sample_table_wedit, grepl("P635L", mutant))
# Transpose count matrices
# t_counts_P635L <- t(featurecounts_P635L$counts)
t_counts_P635L_after_base_editing <- t(featurecounts_P635L_after_base_editing$counts)
# Rownames to columns
# t_counts_P635L <- rownames_to_column(as.data.frame(t_counts_P635L), 
#                                      var = "sample")
t_counts_P635L_after_base_editing <- rownames_to_column(as.data.frame(t_counts_P635L_after_base_editing), 
                                                        var = "sample")
# Add metadata
# t_counts_P635L_wmet <- dplyr::left_join(t_counts_P635L, 
#                                          sample_table_noedit_P635L)
t_counts_P635L_after_base_editing_wmet <- dplyr::left_join(t_counts_P635L_after_base_editing, 
                                                           sample_table_wedit_P635L)
t_counts_P635L_after_base_editing_wmet <- dplyr::filter(t_counts_P635L_after_base_editing_wmet, 
                                                        !is.na(treatment)) # I don't know why but the liver samples are there --> need to check

# Compute mean per group
# mean_t_counts_P635L_wmet <- t_counts_P635L_wmet %>% 
#   group_by(genotype) %>% 
#   summarise(across(contains("ENSMUSG"), mean))
mean_t_counts_P635L_after_base_editing_wmet <- t_counts_P635L_after_base_editing_wmet %>% 
  group_by(treatment) %>% 
  summarise(across(contains("ENSMUSG"), mean))

# Compute correlation
# # Long format
# mean_counts_P635L_wmet <- pivot_longer(mean_t_counts_P635L_wmet, 
#                                        contains("ENSMUSG"), 
#                                        names_to = "gene", 
#                                        values_to = "counts")
# mean_counts_P635L_after_base_editing_wmet <- pivot_longer(mean_t_counts_P635L_after_base_editing_wmet, 
#                                                           contains("ENSMUSG"), 
#                                                           names_to = "gene", 
#                                                           values_to = "counts")
# # Combine
# mean_counts_P635L_wmet$experiment <- "P635L"
# mean_counts_P635L_wmet$treatment <- "none"
# mean_counts_P635L_after_base_editing_wmet$experiment <- "P635L_after_base_editing"
# mean_counts_P635L_combinedexp <- rbind(mean_counts_P635L_wmet, mean_counts_P635L_after_base_editing_wmet)
# #
# mean_counts_P635L_combinedexp %>% 
#   dplyr::filter(genotype == "hom") %>% 
#   group_by(treatment) %>% 
#   summarize()
#   

# Create rownames
t_counts_P635L_after_base_editing_wmet_tomatrix <- unite(t_counts_P635L_after_base_editing_wmet, 
                                                         sample_treat, treatment, run,  sep = "_")
t_counts_P635L_after_base_editing_wmet_tomatrix <- column_to_rownames(t_counts_P635L_after_base_editing_wmet_tomatrix, 
                                                                      "sample_treat")
# Matrix
matrix_t_counts_P635L_after_base_editing_wmet <- dplyr::select(t_counts_P635L_after_base_editing_wmet_tomatrix,
              contains("ENSMUSG"))
matrix_t_counts_P635L_after_base_editing_wmet <- as.matrix(matrix_t_counts_P635L_after_base_editing_wmet)
# Correlation
cor_matrix_t_counts_P635L_after_base_editing_wmet <- 
  cor(t(matrix_t_counts_P635L_after_base_editing_wmet), method = "spearman")
# Pivot longer
df_cor_matrix_t_counts_P635L_after_base_editing_wmet <- rownames_to_column(as.data.frame(cor_matrix_t_counts_P635L_after_base_editing_wmet), 
                                                                           "sample_treat1")
lg_cor_matrix_t_counts_P635L_after_base_editing_wmet <- pivot_longer(df_cor_matrix_t_counts_P635L_after_base_editing_wmet, 
             cols = !sample_treat1,
             names_to = "sample_treat2")
             
# Plot
plot_cor_P635L_after_base_editing_wmet <- ggplot(data = lg_cor_matrix_t_counts_P635L_after_base_editing_wmet, 
       aes(x=sample_treat1, y=sample_treat2, fill=value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_continuous(high = "#132B43", low = "#56B1F7") +
  ggtitle("Spearman correlation P635L after base editing")
  # scale_fill_continuous(trans = 'reverse')
ggsave("results/plot_cor_P635L_after_base_editing_wmet.pdf", 
       plot_cor_P635L_after_base_editing_wmet,
       height = 10, width = 10)

# Exclude some samples
# PBS_22s002888, Nterm_SpRY_and_Cterm_SpRY_gRNA5_22s002878 and  Nterm_SpRY_and_Cterm_SpRY_gRNA5_22s002879
t_counts_P635L_after_base_editing_wmet_tomatrix_norownames <- unite(t_counts_P635L_after_base_editing_wmet, 
                                                         sample_treat, treatment, run,  sep = "_")
t_counts_P635L_after_base_editing_wmet_tomatrix_norownames_filtered <- dplyr::filter(t_counts_P635L_after_base_editing_wmet_tomatrix_norownames,
              !(sample_treat %in% c("PBS_22s002888", "Nterm_SpRY_and_Cterm_SpRY_gRNA5_22s002878", "Nterm_SpRY_and_Cterm_SpRY_gRNA5_22s002879")))
t_counts_P635L_after_base_editing_wmet_tomatrix_filtered <- column_to_rownames(t_counts_P635L_after_base_editing_wmet_tomatrix_norownames_filtered, 
                                                                      "sample_treat")
# Matrix
matrix_t_counts_P635L_after_base_editing_wmet_filtered <- dplyr::select(t_counts_P635L_after_base_editing_wmet_tomatrix_filtered,
                                                               contains("ENSMUSG"))
matrix_t_counts_P635L_after_base_editing_wmet_filtered <- as.matrix(matrix_t_counts_P635L_after_base_editing_wmet_filtered)
# Correlation
cor_matrix_t_counts_P635L_after_base_editing_wmet_filtered <- 
  cor(t(matrix_t_counts_P635L_after_base_editing_wmet_filtered), method = "spearman")
# Pivot longer
df_cor_matrix_t_counts_P635L_after_base_editing_wmet_filtered <- rownames_to_column(as.data.frame(cor_matrix_t_counts_P635L_after_base_editing_wmet_filtered), 
                                                                           "sample_treat1")
lg_cor_matrix_t_counts_P635L_after_base_editing_wmet_filtered <- pivot_longer(df_cor_matrix_t_counts_P635L_after_base_editing_wmet_filtered, 
                                                                     cols = !sample_treat1,
                                                                     names_to = "sample_treat2")

# Plot
plot_cor_P635L_after_base_editing_wmet_filtered <- ggplot(data = lg_cor_matrix_t_counts_P635L_after_base_editing_wmet_filtered, 
                                                 aes(x=sample_treat1, y=sample_treat2, fill=value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_continuous(high = "#132B43", low = "#56B1F7") +
  ggtitle("Spearman correlation P635L after base editing - 3 samples filtered")
# scale_fill_continuous(trans = 'reverse')
ggsave("results/plot_cor_P635L_after_base_editing_wmet_filtered.pdf", 
       plot_cor_P635L_after_base_editing_wmet_filtered,
       height = 10, width = 10)


# R636Q
# Load the data
# After base editing
featurecounts_R636Q_after_base_editing <- 
  readRDS("results/robject_featurecounts_R636Q_after_base_editing.rds")

# Add metadata
# Sample table
sample_table_wedit <- read.table("/g/steinmetz/project/dcm_lgreads/mouse_bulkRNA/config/samples_afterbaseediting_metadata.txt",
                                 header = TRUE)
# Add info
# Simplify sample table
sample_table_wedit$sample <- paste("results/STAR", 
                                   sample_table_wedit$mutant, 
                                   sample_table_wedit$run, 
                                   "Aligned.sortedByCoord.out.bam", sep = "/")
# Filter
sample_table_wedit_R636Q <- dplyr::filter(sample_table_wedit, grepl("R636Q", mutant))
# Transpose count matrices
t_counts_R636Q_after_base_editing <- t(featurecounts_R636Q_after_base_editing$counts)
# Rownames to columns
t_counts_R636Q_after_base_editing <- rownames_to_column(as.data.frame(t_counts_R636Q_after_base_editing), 
                                                        var = "sample")
# Add metadata
t_counts_R636Q_after_base_editing_wmet <- dplyr::left_join(t_counts_R636Q_after_base_editing, 
                                                           sample_table_wedit_R636Q)
t_counts_R636Q_after_base_editing_wmet <- dplyr::filter(t_counts_R636Q_after_base_editing_wmet, 
                                                        !is.na(treatment)) # I don't know why but the liver samples are there --> need to check

# # Compute mean per group
# mean_t_counts_R636Q_after_base_editing_wmet <- t_counts_R636Q_after_base_editing_wmet %>% 
#   group_by(treatment) %>% 
#   summarise(across(contains("ENSMUSG"), mean))

# Compute correlation
# Create rownames
t_counts_R636Q_after_base_editing_wmet_tomatrix <- unite(t_counts_R636Q_after_base_editing_wmet, 
                                                         sample_treat, treatment, run,  sep = "_")
t_counts_R636Q_after_base_editing_wmet_tomatrix <- column_to_rownames(t_counts_R636Q_after_base_editing_wmet_tomatrix, 
                                                                      "sample_treat")
# Matrix
matrix_t_counts_R636Q_after_base_editing_wmet <- dplyr::select(t_counts_R636Q_after_base_editing_wmet_tomatrix,
                                                               contains("ENSMUSG"))
matrix_t_counts_R636Q_after_base_editing_wmet <- as.matrix(matrix_t_counts_R636Q_after_base_editing_wmet)
# Correlation
cor_matrix_t_counts_R636Q_after_base_editing_wmet <- 
  cor(t(matrix_t_counts_R636Q_after_base_editing_wmet), method = "spearman")
# Pivot longer
df_cor_matrix_t_counts_R636Q_after_base_editing_wmet <- rownames_to_column(as.data.frame(cor_matrix_t_counts_R636Q_after_base_editing_wmet), 
                                                                           "sample_treat1")
lg_cor_matrix_t_counts_R636Q_after_base_editing_wmet <- pivot_longer(df_cor_matrix_t_counts_R636Q_after_base_editing_wmet, 
                                                                     cols = !sample_treat1,
                                                                     names_to = "sample_treat2")

# Plot
plot_cor_R636Q_after_base_editing_wmet <- ggplot(data = lg_cor_matrix_t_counts_R636Q_after_base_editing_wmet, 
                                                 aes(x=sample_treat1, y=sample_treat2, fill=value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_continuous(high = "#132B43", low = "#56B1F7") +
  ggtitle("Spearman correlation R636Q after base editing")
# scale_fill_continuous(trans = 'reverse')
ggsave("results/plot_cor_R636Q_after_base_editing_wmet.pdf", 
       plot_cor_R636Q_after_base_editing_wmet,
       height = 10, width = 10)
