# DESeq2 analysis with STAR alignment

# Load libraries
library(DESeq2)
library(Rsubread)
library(tidyverse)
library(ggplot2)
library(ggfortify)
library(pheatmap)
library(rtracklayer)

# Read gtf file
# gen_vM29_annot <- read.delim("resources/gencode.vM29.annotation.gtf", header=F, comment.char="#")
gen_vM29_annot <- import("resources/gencode.vM29.annotation.gtf")
# Extract info
gen_vM29_annot_gene <- gen_vM29_annot[,c("gene_id", "gene_name")]
gen_vM29_annot_gene_uniq <- unique(as.data.frame(gen_vM29_annot_gene)[,c("gene_id", "gene_name")])
colnames(gen_vM29_annot_gene_uniq)[1] <- c("gene")

# Load count matrices
load("results/robject_featurecounts_allsamples.rda")
# featurecounts_allsamples
# Sample table
sample_table <- read.table("/g/steinmetz/project/dcm_lgreads/mouse_bulkRNA/config/samples_metadata.txt",
                           header = TRUE)
# read.table(snakemake@input$meta, header = TRUE)
# Simplify sample table
sample_table$sample <- paste("results/STAR", 
                             sample_table$mutant, 
                             sample_table$run, 
                             "Aligned.sortedByCoord.out.bam", sep = "/")

# PCA
# Prepare data for PCA
counts_all_mice <- t(featurecounts_allsamples$counts)
counts_all_mice_wnames <- rownames_to_column(as.data.frame(counts_all_mice), var = "sample")
counts_all_mice_wmeta <- left_join(counts_all_mice_wnames, sample_table)
# Perform PCA
# All mice
pca_all_mice <- prcomp(counts_all_mice)
plot_pca_all_mice <- autoplot(pca_all_mice,
                              data = counts_all_mice_wmeta,
                              colour = "genotype",
                              shape = "mutant", 
                              size = 3,
                              alpha = 0.8, 
                              label = TRUE,  
                              label.label = "run",
                              label.position = position_nudge(x = 0, y = 0.015),
                              label.size = 2) +
  theme_bw() +
  ggtitle("PCA with all mice and all genes")
ggsave("results/plot_STARalignment_pca_all_mice_genes.pdf",
       plot = plot_pca_all_mice,
       width = 7, height = 5, units = "in")

# DESeq2
# DESeqDataSetFromMatrix
# P635L
# Filter 
sample_table_P635L <- dplyr::filter(sample_table, mutant == "P635L")
# Extract counts
featurecounts_P635L <- featurecounts_allsamples$counts[,sample_table_P635L$sample]
# Check order
if(all(colnames(featurecounts_P635L) == sample_table_P635L$sample)){
  colnames(featurecounts_P635L) <- sample_table_P635L$run
} else {
  # Rename columns
  stop(paste0("Check that the order of the columns of the counts for P635L mutant mice corresponds to the sample table"))
}
# Create deseq2 object
dds_P635L <- DESeqDataSetFromMatrix(countData = featurecounts_P635L,
                                    colData = sample_table_P635L[,c("run", "genotype")],
                                    design = ~ genotype)
# Perform DE test
dds_P635L <- DESeq(dds_P635L)
# P635L hom vs WT
res_P635L_homvswt <- results(dds_P635L, contrast=c("genotype","hom","wt"))
res_P635L_homvswt_sig <- dplyr::filter(as.data.frame(res_P635L_homvswt), 
                                       padj < 0.05)
res_P635L_homvswt_sig <- rownames_to_column(res_P635L_homvswt_sig, var = "gene")
write.table(x = res_P635L_homvswt_sig,
            file = "results/table_STARalignment_res_P635L_homvswt_sig.txt",
            row.names = FALSE,
            col.names = TRUE)

# R636Q
# Filter 
sample_table_R636Q <- dplyr::filter(sample_table, mutant == "R636Q")
# Extract counts
featurecounts_R636Q <- featurecounts_allsamples$counts[,sample_table_R636Q$sample]
# Check order
if(all(colnames(featurecounts_R636Q) == sample_table_R636Q$sample)){
  colnames(featurecounts_R636Q) <- sample_table_R636Q$run
} else {
  # Rename columns
  stop(paste0("Check that the order of the columns of the counts for R636Q mutant mice corresponds to the sample table"))
}
# Create deseq2 object
dds_R636Q <- DESeqDataSetFromMatrix(countData = featurecounts_R636Q,
                                    colData = sample_table_R636Q[,c("run", "genotype")],
                                    design = ~ genotype)
# Perform DE test
dds_R636Q <- DESeq(dds_R636Q)
# R636Q hom vs WT
res_R636Q_homvswt <- results(dds_R636Q, contrast=c("genotype","hom","wt"))
res_R636Q_homvswt_sig <- dplyr::filter(as.data.frame(res_R636Q_homvswt), 
                                       padj < 0.05)
res_R636Q_homvswt_sig <- rownames_to_column(res_R636Q_homvswt_sig, var = "gene")
write.table(x = res_R636Q_homvswt_sig,
            file = "results/table_STARalignment_res_R636Q_homvswt_sig.txt",
            row.names = FALSE,
            col.names = TRUE)

# P635L vs R636Q wt
# Filter 
sample_table_wt <- dplyr::filter(sample_table, genotype == "wt")
# Extract counts
featurecounts_wt <- featurecounts_allsamples$counts[,sample_table_wt$sample]
# Check order
if(all(colnames(featurecounts_wt) == sample_table_wt$sample)){
  colnames(featurecounts_wt) <- sample_table_wt$run
} else {
  # Rename columns
  stop(paste0("Check that the order of the columns of the counts for WT mice corresponds to the sample table"))
}
# Create deseq2 object
dds_wt <- DESeqDataSetFromMatrix(countData = featurecounts_wt,
                                 colData = sample_table_wt[,c("mutant", "run")],
                                 design = ~ mutant)
# Perform DE test
dds_wt <- DESeq(dds_wt)
# P635L vs R636Q wt
res_P635LvsR636Q_wt <- results(dds_wt, contrast=c("mutant", "P635L", "R636Q"))
res_P635LvsR636Q_wt_sig <- dplyr::filter(as.data.frame(res_P635LvsR636Q_wt), 
                                         padj < 0.05)
res_P635LvsR636Q_wt_sig <- rownames_to_column(res_P635LvsR636Q_wt_sig, var = "gene")
write.table(x = res_P635LvsR636Q_wt_sig,
            file = "results/table_STARalignment_res_P635LvsR636Q_wt_sig.txt",
            row.names = FALSE,
            col.names = TRUE)

# Heatmap function
# mutation: character, mutation
# dds: DESeqDataSet, DESeq2 object
# res: data.frame, significant results from DESeq DE test
# sample_tab: data.frame, table with genotype for each sample
# annot: data.frame, one column contains the ensembl IDs (called gene) and the other the gene_name
# filter_col: character, column to filter on
# gen: vector, selected phenotype
# top: numeric, top gene to display
function_my_heatmap <- function(mutation, dds, res, sample_tab, annot, filter_col, gen, top){
  # Order genes based on adjusted p-values
  res_ordered <- res[order(res$padj),]
  # Select genotype samples
  df_table <- dplyr::filter(sample_tab, 
                            get(filter_col) %in% gen) %>% 
    dplyr::select(run, filter_col)
  row.names(df_table) <- df_table$run
  # Extract top genes
  counts_allsamples <- counts(dds, normalized=TRUE)[res_ordered[1:top, "gene"],]
  counts_allsamples <- rownames_to_column(as.data.frame(counts_allsamples), 
                                          var = "gene")
  # Select genotype samples
  counts_sig <- dplyr::select(counts_allsamples, c(df_table$run, "gene"))
  # Add gene name
  counts_sig_wgname <- dplyr::left_join(counts_sig, annot)
  counts_sig_wgname <- column_to_rownames(counts_sig_wgname, "gene_name")
  # Plot
  pheatmap(counts_sig_wgname[,1:length(df_table$run)], 
           cluster_rows = TRUE, 
           show_rownames = TRUE,
           show_colnames = FALSE,
           cluster_cols = TRUE,
           annotation_col = df_table,
           scale = "row",
           main = mutation)
}
pdf("results/heatmap_STARalignment_R636Q_homvswt_top20.pdf")
function_my_heatmap(mutation = "R636Q",
                    dds = dds_R636Q,
                    res = res_R636Q_homvswt_sig,
                    sample_tab = sample_table_R636Q, 
                    annot = gen_vM29_annot_gene_uniq,
                    gen = c("hom", "wt"),
                    filter_col = "genotype",
                    top = 20)
dev.off()
pdf("results/heatmap_STARalignment_P635L_homvswt_top20.pdf")
function_my_heatmap(mutation = "P635L",
                    dds = dds_P635L,
                    res = res_P635L_homvswt_sig,
                    sample_tab = sample_table_P635L, 
                    annot = gen_vM29_annot_gene_uniq,
                    gen = c("hom", "wt"),
                    filter_col = "genotype",
                    top = 20)
dev.off()
pdf("results/heatmap_STARalignment_P635LvsR636Q_wt_top20.pdf")
function_my_heatmap(mutation = "WT",
                    dds = dds_wt,
                    res = res_P635LvsR636Q_wt_sig,
                    sample_tab = sample_table_wt, 
                    annot = gen_vM29_annot_gene_uniq,
                    gen = c("P635L", "R636Q"),
                    filter_col = "mutant",
                    top = 20)
dev.off()

# Compute fold change per individuals for each mutation
# log2(hom/wt)
# P635L
dim(counts(dds_P635L, normalized = TRUE))
geneexp_dds_P635L <- t(counts(dds_P635L, normalized = TRUE))
# Average per gene for WT
wt_IDs <- dplyr::filter(sample_table_P635L, genotype == "wt") %>% pull(run)
pergene_P635L_WT <- dplyr::select(as.data.frame(counts(dds_P635L, , normalized = TRUE)), wt_IDs)
avg_pergene_P635L_WT <- rowMeans(pergene_P635L_WT)
# Remove WT columns
geneexp_dds_P635L_hom_het <- dplyr::select(as.data.frame(counts(dds_P635L, normalized = TRUE)), -wt_IDs)
# Long format
geneexp_dds_P635L_hom_het <- rownames_to_column(.data = geneexp_dds_P635L_hom_het, 
                                                var = "gene")
geneexp_dds_P635L_hom_het_lg <- pivot_longer(geneexp_dds_P635L_hom_het, 
                                             cols = !gene, 
                                             names_to = "run")
# Join genotype
geneexp_dds_P635L_hom_het_lg_wgenot <- left_join(geneexp_dds_P635L_hom_het_lg, 
                                                 sample_table_P635L[, c("mutant", "run", "genotype")])
# Join wt average value
avg_pergene_P635L_WT_df <- data.frame(gene=names(avg_pergene_P635L_WT), 
                                      wt_avg=avg_pergene_P635L_WT, 
                                      row.names=NULL)
geneexp_dds_P635L_hom_het_lg_wgenot_wwt <- left_join(geneexp_dds_P635L_hom_het_lg_wgenot, 
                                                     avg_pergene_P635L_WT_df)
# Filter genes with zero reads in WT mice as I won't be able to compute a log2 FC for these genes
geneexp_dds_P635L_hom_het_lg_wgenot_wwt <- dplyr::filter(geneexp_dds_P635L_hom_het_lg_wgenot_wwt,
                                                         wt_avg != 0)
# FC
geneexp_dds_P635L_hom_het_lg_log2FC <- geneexp_dds_P635L_hom_het_lg_wgenot_wwt %>% 
  group_by(gene, run) %>% 
  mutate(log2_FC = log2(value/wt_avg))
# Replace -Inf by NA
# dplyr::filter(geneexp_dds_P635L_hom_het_lg_log2FC, log2_FC == -Inf)
geneexp_dds_P635L_hom_het_lg_log2FC <- mutate(geneexp_dds_P635L_hom_het_lg_log2FC, 
                                              log2_FC = replace(log2_FC, log2_FC == -Inf, NA))

# R636Q
dim(counts(dds_R636Q, normalized = TRUE))
geneexp_dds_R636Q <- t(counts(dds_R636Q, normalized = TRUE))
# Average per gene for WT
wt_IDs <- dplyr::filter(sample_table_R636Q, genotype == "wt") %>% pull(run)
pergene_R636Q_WT <- dplyr::select(as.data.frame(counts(dds_R636Q, normalized = TRUE)), wt_IDs)
avg_pergene_R636Q_WT <- rowMeans(pergene_R636Q_WT)
# Remove WT columns
geneexp_dds_R636Q_hom_het <- dplyr::select(as.data.frame(counts(dds_R636Q, normalized = TRUE)), -wt_IDs)
# Long format
geneexp_dds_R636Q_hom_het <- rownames_to_column(.data = geneexp_dds_R636Q_hom_het, 
                                                var = "gene")
geneexp_dds_R636Q_hom_het_lg <- pivot_longer(geneexp_dds_R636Q_hom_het, 
                                             cols = !gene, 
                                             names_to = "run")
# Join genotype
geneexp_dds_R636Q_hom_het_lg_wgenot <- left_join(geneexp_dds_R636Q_hom_het_lg, 
                                                 sample_table_R636Q[, c("mutant", "run", "genotype")])
# Join wt average value
avg_pergene_R636Q_WT_df <- data.frame(gene=names(avg_pergene_R636Q_WT), 
                                      wt_avg=avg_pergene_R636Q_WT, 
                                      row.names=NULL)
geneexp_dds_R636Q_hom_het_lg_wgenot_wwt <- left_join(geneexp_dds_R636Q_hom_het_lg_wgenot, 
                                                     avg_pergene_R636Q_WT_df)
# Filter genes with zero reads in WT mice as I won't be able to compute a log2 FC for these genes
geneexp_dds_R636Q_hom_het_lg_wgenot_wwt <- dplyr::filter(geneexp_dds_R636Q_hom_het_lg_wgenot_wwt,
                                                         wt_avg != 0)
# FC
geneexp_dds_R636Q_hom_het_lg_log2FC <- geneexp_dds_R636Q_hom_het_lg_wgenot_wwt %>% 
  group_by(gene, run) %>% 
  mutate(log2_FC = log2(value/wt_avg))
# Replace -Inf by NA
# dplyr::filter(geneexp_dds_R636Q_hom_het_lg_log2FC, log2_FC == -Inf)
geneexp_dds_R636Q_hom_het_lg_log2FC <- mutate(geneexp_dds_R636Q_hom_het_lg_log2FC, 
                                              log2_FC = replace(log2_FC, log2_FC == -Inf, NA))

# Heatmap
# Filter gene list based on Markus list
DEgenelist_p_inf_1e4 <- read.table("results/DEgenelist_p_inf_1e-4.txt")
DEgenelist_p_inf_1e6 <- read.table("results/DEgenelist_p_inf_1e-6.txt")
# Merge with ensembl IDs
DEgenelist_p_inf_1e4_wname <- left_join(data.frame("gene_name" = DEgenelist_p_inf_1e4[,1]), 
                                        gen_vM29_annot_gene_uniq)
DEgenelist_p_inf_1e6_wname <- left_join(data.frame("gene_name" = DEgenelist_p_inf_1e6[,1]), 
                                        gen_vM29_annot_gene_uniq)
# Merge gene name with fold change and order genes according to the avg log2 FC
fct_merge_FC_wgenename_and_order <- function(geneexp, DEgenelist){
  # Merge gene name with fold change
  df_merge <- left_join(geneexp, DEgenelist)
  # Filter NA row 
  df_filtered <- dplyr::filter(df_merge, 
                               run != "NA")
  # Order per log2FC average
  df_filtered_avglg2fc <- df_filtered %>% 
    group_by(gene) %>% 
    mutate(avg_log2FC = mean(log2_FC, na.rm = TRUE))
  df_filtered_avglg2fc[order(df_filtered_avglg2fc$avg_log2FC, decreasing = TRUE),]
}
# P635L p < 1e-4
heatmap_data_P635L_hom_het_p_inf_1e4 <- fct_merge_FC_wgenename_and_order(
  geneexp = DEgenelist_p_inf_1e4_wname,
  DEgenelist = geneexp_dds_P635L_hom_het_lg_log2FC
)
# P635L p < 1e-6
heatmap_data_P635L_hom_het_p_inf_1e6 <- fct_merge_FC_wgenename_and_order(
  geneexp = DEgenelist_p_inf_1e6_wname,
  DEgenelist = geneexp_dds_P635L_hom_het_lg_log2FC
)
# R636Q p < 1e-4
heatmap_data_R636Q_hom_het_p_inf_1e4 <- fct_merge_FC_wgenename_and_order(
  geneexp = DEgenelist_p_inf_1e4_wname,
  DEgenelist = geneexp_dds_R636Q_hom_het_lg_log2FC
)
# R636Q p < 1e-6
heatmap_data_R636Q_hom_het_p_inf_1e6 <- fct_merge_FC_wgenename_and_order(
  geneexp = DEgenelist_p_inf_1e6_wname,
  DEgenelist = geneexp_dds_R636Q_hom_het_lg_log2FC
)

# Transform to wide format + Collapse run with genotype
fct_wide_gen_order <- function(heatmapdata){
  # Wide format
  df_w <- pivot_wider(heatmapdata, 
                      id_cols = c("run", "genotype"),
                      names_from = "gene_name", 
                      values_from = "log2_FC")
  # Collapse run with genotype
  df_w_gen <- df_w %>% 
    rowwise() %>% 
    mutate(sample = paste(run, genotype, sep = "_"), .before = run)
  df_w_gen <- dplyr::select(df_w_gen, -c(run, genotype))
  column_to_rownames(df_w_gen, var = "sample")
}
# P635L p < 1e-4
heatmap_data_P635L_hom_het_p_inf_1e4_w <- fct_wide_gen_order(
  heatmapdata = heatmap_data_P635L_hom_het_p_inf_1e4
)
# P635L p < 1e-6
heatmap_data_P635L_hom_het_p_inf_1e6_w <- fct_wide_gen_order(
  heatmapdata = heatmap_data_P635L_hom_het_p_inf_1e6
)
# R636Q p < 1e-4
heatmap_data_R636Q_hom_het_p_inf_1e4_w <- fct_wide_gen_order(
  heatmapdata = heatmap_data_R636Q_hom_het_p_inf_1e4
)
# R636Q p < 1e-6
heatmap_data_R636Q_hom_het_p_inf_1e6_w <- fct_wide_gen_order(
  heatmapdata = heatmap_data_R636Q_hom_het_p_inf_1e6
)

# Heatmaps
png("results/heatmap_STARalignment_P635L_hom_het_p_inf_1e4.png", 
    width = 8, height = 20, units = "in", res = 1000)
pheatmap(t(heatmap_data_P635L_hom_het_p_inf_1e4_w),
         na_col = "black",
         cluster_rows = FALSE, 
         show_rownames = TRUE,
         show_colnames = TRUE,
         cluster_cols = FALSE,
         fontsize_row = 2, 
         main = "P635L with DE genes p < 1e-4")
dev.off()
png("results/heatmap_STARalignment_P635L_hom_het_p_inf_1e6.png", 
    width = 8, height = 20, units = "in", res = 1000)
pheatmap(t(heatmap_data_P635L_hom_het_p_inf_1e6_w),
         na_col = "black",
         cluster_rows = FALSE, 
         show_rownames = TRUE,
         show_colnames = TRUE,
         cluster_cols = FALSE,
         fontsize_row = 2, 
         main = "P635L with DE genes p < 1e-6")
dev.off()
# Order heatmap rows
idx_het <- grep("*het$", rownames(heatmap_data_R636Q_hom_het_p_inf_1e4_w))
idx_hom <- grep("*hom$", rownames(heatmap_data_R636Q_hom_het_p_inf_1e4_w))
idx <- c(idx_hom, idx_het)
png("results/heatmap_STARalignment_R636Q_hom_het_p_inf_1e4.png", 
    width = 8, height = 20, units = "in", res = 1000)
pheatmap(t(heatmap_data_R636Q_hom_het_p_inf_1e4_w[idx,]),
         na_col = "black",
         cluster_rows = FALSE, 
         show_rownames = TRUE,
         show_colnames = TRUE,
         cluster_cols = FALSE,
         fontsize_row = 2, 
         main = "R636Q with DE genes p < 1e-4")
dev.off()
png("results/heatmap_STARalignment_R636Q_hom_het_p_inf_1e6.png", 
    width = 8, height = 20, units = "in", res = 1000)
pheatmap(t(heatmap_data_R636Q_hom_het_p_inf_1e6_w[idx,]),
         na_col = "grey90",
         cluster_rows = FALSE, 
         show_rownames = TRUE,
         show_colnames = TRUE,
         cluster_cols = FALSE,
         fontsize_row = 2, 
         main = "R636Q with DE genes p < 1e-6")
dev.off()
