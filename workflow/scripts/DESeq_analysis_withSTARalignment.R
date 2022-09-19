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
# P635L HET vs WT
res_P635L_hetvswt <- results(dds_P635L, contrast=c("genotype","het","wt"))
res_P635L_hetvswt_sig <- dplyr::filter(as.data.frame(res_P635L_hetvswt), 
                                       padj < 0.05)
res_P635L_hetvswt_sig <- rownames_to_column(res_P635L_hetvswt_sig, var = "gene")
write.table(x = res_P635L_hetvswt_sig,
            file = "results/table_STARalignment_res_P635L_hetvswt_sig.txt",
            row.names = FALSE,
            col.names = TRUE)
res_P635L_hetvswt_sig_withgenename <- dplyr::left_join(res_P635L_hetvswt_sig, gen_vM29_annot_gene_uniq)
write.table(x = res_P635L_hetvswt_sig_withgenename,
            file = "results/table_STARalignment_res_P635L_hetvswt_sig_withgenename.txt",
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
# R636Q het vs WT
res_R636Q_hetvswt <- results(dds_R636Q, contrast=c("genotype","het","wt"))
res_R636Q_hetvswt_sig <- dplyr::filter(as.data.frame(res_R636Q_hetvswt), 
                                       padj < 0.05)
res_R636Q_hetvswt_sig <- rownames_to_column(res_R636Q_hetvswt_sig, var = "gene")
write.table(x = res_R636Q_hetvswt_sig,
            file = "results/table_STARalignment_res_R636Q_hetvswt_sig.txt",
            row.names = FALSE,
            col.names = TRUE)
res_R636Q_hetvswt_sig_withgenename <- dplyr::left_join(res_R636Q_hetvswt_sig, gen_vM29_annot_gene_uniq)
write.table(x = res_R636Q_hetvswt_sig_withgenename,
            file = "results/table_STARalignment_res_R636Q_hetvswt_sig_withgenename.txt",
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
pdf("results/heatmap_STARalignment_R636Q_hetvswt_top20.pdf")
function_my_heatmap(mutation = "R636Q",
                    dds = dds_R636Q,
                    res = res_R636Q_hetvswt_sig,
                    sample_tab = sample_table_R636Q, 
                    annot = gen_vM29_annot_gene_uniq,
                    gen = c("het", "wt"),
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
pdf("results/heatmap_STARalignment_P635L_hetvswt_top20.pdf")
function_my_heatmap(mutation = "P635L",
                    dds = dds_P635L,
                    res = res_P635L_hetvswt_sig,
                    sample_tab = sample_table_P635L, 
                    annot = gen_vM29_annot_gene_uniq,
                    gen = c("het", "wt"),
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
pergene_P635L_WT <- dplyr::select(as.data.frame(counts(dds_P635L, normalized = TRUE)), wt_IDs)
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
# Extract list of genes with zero reads in WT and combine with gene name
geneexp_dds_P635L_hom_het_lg_wgenot_wwt_zeroreads <- dplyr::filter(geneexp_dds_P635L_hom_het_lg_wgenot_wwt,
                                                                   wt_avg == 0) %>% 
  pull(gene) %>% 
  unique
geneexp_dds_P635L_hom_het_lg_wgenot_wwt_zeroreads <- left_join(data.frame(gene = geneexp_dds_P635L_hom_het_lg_wgenot_wwt_zeroreads), 
                                                               gen_vM29_annot_gene_uniq)
write.table(geneexp_dds_P635L_hom_het_lg_wgenot_wwt_zeroreads,
            file = "results/geneexp_P635L_hom_het_WT_zeroreads.txt", 
            row.names = FALSE,
            col.names = TRUE)
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
# Extract list of genes with zero reads in WT
geneexp_dds_R636Q_hom_het_lg_wgenot_wwt_zeroreads <- dplyr::filter(geneexp_dds_R636Q_hom_het_lg_wgenot_wwt,
                                                                   wt_avg == 0) %>% 
  pull(gene) %>% 
  unique
geneexp_dds_R636Q_hom_het_lg_wgenot_wwt_zeroreads <- left_join(data.frame(gene = geneexp_dds_R636Q_hom_het_lg_wgenot_wwt_zeroreads), 
                                                               gen_vM29_annot_gene_uniq)
write.table(geneexp_dds_R636Q_hom_het_lg_wgenot_wwt_zeroreads,
            file = "results/geneexp_R636Q_hom_het_WT_zeroreads.txt", 
            row.names = FALSE,
            col.names = TRUE)
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

# Combine WT gene lists with zero genes
geneexp_dds_P635LR636Q_hom_het_lg_wgenot_wwt_zeroreads_combined <- 
  rbind(data.frame(geneexp_dds_R636Q_hom_het_lg_wgenot_wwt_zeroreads, "mutant" = "R636Q"),
        data.frame(geneexp_dds_P635L_hom_het_lg_wgenot_wwt_zeroreads, "mutant" = "P635L"))
geneexp_dds_P635LR636Q_hom_het_lg_wgenot_wwt_zeroreads_combined_short <- 
  geneexp_dds_P635LR636Q_hom_het_lg_wgenot_wwt_zeroreads_combined %>% 
  group_by(gene_name) %>% 
  summarize(mutants = paste0(unique(mutant), collapse = "_"))
write.table(geneexp_dds_P635LR636Q_hom_het_lg_wgenot_wwt_zeroreads_combined_short,
            file = "results/geneexp_P635LR636Q_hom_het_WT_zeroreads.txt", 
            row.names = FALSE,
            col.names = TRUE)

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
idx_R636Q_het <- grep("*het$", rownames(heatmap_data_R636Q_hom_het_p_inf_1e4_w))
idx_R636Q_hom <- grep("*hom$", rownames(heatmap_data_R636Q_hom_het_p_inf_1e4_w))
idx_R636Q <- c(idx_R636Q_hom, idx_R636Q_het)
png("results/heatmap_STARalignment_R636Q_hom_het_p_inf_1e4.png", 
    width = 8, height = 20, units = "in", res = 1000)
pheatmap(t(heatmap_data_R636Q_hom_het_p_inf_1e4_w[idx_R636Q,]),
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
pheatmap(t(heatmap_data_R636Q_hom_het_p_inf_1e6_w[idx_R636Q,]),
         na_col = "grey90",
         cluster_rows = FALSE, 
         show_rownames = TRUE,
         show_colnames = TRUE,
         cluster_cols = FALSE,
         fontsize_row = 2, 
         main = "R636Q with DE genes p < 1e-6")
dev.off()

# All individuals on one heatmap
# p < 1e-4
heatmap_data_P635LR636Q_hom_het_p_inf_1e4_w <- 
  merge(heatmap_data_P635L_hom_het_p_inf_1e4_w, 
        heatmap_data_R636Q_hom_het_p_inf_1e4_w[idx_R636Q,],
        all = TRUE, 
        sort = FALSE)
rownames(heatmap_data_P635LR636Q_hom_het_p_inf_1e4_w) <-
  c(rownames(heatmap_data_P635L_hom_het_p_inf_1e4_w),
    rownames(heatmap_data_R636Q_hom_het_p_inf_1e4_w[idx_R636Q,]))
# p < 1e-6
heatmap_data_P635LR636Q_hom_het_p_inf_1e6_w <- 
  merge(heatmap_data_P635L_hom_het_p_inf_1e6_w, 
        heatmap_data_R636Q_hom_het_p_inf_1e6_w[idx_R636Q,],
        all = TRUE, 
        sort = FALSE)
rownames(heatmap_data_P635LR636Q_hom_het_p_inf_1e6_w) <-
  c(rownames(heatmap_data_P635L_hom_het_p_inf_1e6_w),
    rownames(heatmap_data_R636Q_hom_het_p_inf_1e6_w[idx_R636Q,]))
# Annotation
annot_P635LR636Q <- data.frame(
  "genotype" = rep(rep(c("homozygote", "heterozygote"), each = 5), 2), 
  "mutation" = rep(c("P635L", "R636Q"), each = 10)
)
rownames(annot_P635LR636Q) <- rownames(heatmap_data_P635LR636Q_hom_het_p_inf_1e4_w)
png("results/heatmap_STARalignment_P635LR636Q_hom_het_p_inf_1e4.png", 
    width = 8, height = 20, units = "in", res = 1000)
pheatmap(t(heatmap_data_P635LR636Q_hom_het_p_inf_1e4_w),
         na_col = "black",
         cluster_rows = FALSE, 
         show_rownames = TRUE,
         show_colnames = TRUE,
         cluster_cols = FALSE,
         fontsize_row = 2, 
         main = "R636Q and P635L with DE genes p < 1e-4",
         annotation_col = annot_P635LR636Q)
dev.off()
png("results/heatmap_STARalignment_P635LR636Q_hom_het_p_inf_1e6.png", 
    width = 8, height = 20, units = "in", res = 1000)
pheatmap(t(heatmap_data_P635LR636Q_hom_het_p_inf_1e6_w),
         na_col = "black",
         cluster_rows = FALSE, 
         show_rownames = TRUE,
         show_colnames = TRUE,
         cluster_cols = FALSE,
         fontsize_row = 2, 
         main = "R636Q and P635L with DE genes p < 1e-6",
         annotation_col = annot_P635LR636Q)
dev.off()

# Dotplot
genes_of_interest <- 
  c("Nppa", "Casq1", "Tnni2", "Mybpc2", "Tnnt3", 
    "Strit1", "Myh7b", "Hopx", "Aqp4", "Nppb")
# 
# heatmap_data_P635L_gene_int <- 
#   dplyr::filter(heatmap_data_P635L_hom_het_p_inf_1e6,
#               gene_name %in% genes_of_interest)
# heatmap_data_R636Q_gene_int <- 
#   dplyr::filter(heatmap_data_R636Q_hom_het_p_inf_1e6,
#                 gene_name %in% genes_of_interest)
# 
# ggplot(rbind(heatmap_data_P635L_gene_int, heatmap_data_R636Q_gene_int), 
#        aes(x = gene_name, y = run)) +
#   geom_point(aes(size = log1p(value), fill = log2_FC), pch = 21) +
#   # Size with outline
#   # blue to red
#   # geom_point(aes(size = log2_FC, colour = log1p(value))) +
#   scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white",
#                          high = "red", space ="Lab" ) +
#   theme_bw()
# 

df_res_P635L_homvswt <- rownames_to_column(as.data.frame(res_P635L_homvswt), "gene")
full_res_P635L_homvswt_log2FC <- inner_join(df_res_P635L_homvswt, 
                                            heatmap_data_P635L_gene_int)
df_res_R636Q_homvswt <- rownames_to_column(as.data.frame(res_R636Q_homvswt), "gene")
full_res_R636Q_homvswt_log2FC <- inner_join(df_res_R636Q_homvswt, 
                                            heatmap_data_R636Q_gene_int)
# Combine
full_res_P635LR636Q_homvswt_log2FC <- 
  rbind(full_res_P635L_homvswt_log2FC, full_res_R636Q_homvswt_log2FC)
# Add significan column
full_res_P635LR636Q_homvswt_log2FC_sig <-
  full_res_P635LR636Q_homvswt_log2FC %>% 
  mutate(significant = padj < 1e-6)
# Combine mutant and genotype
full_res_P635LR636Q_homvswt_log2FC_sig <- 
  full_res_P635LR636Q_homvswt_log2FC_sig %>% 
  mutate(mutant_gen = paste(mutant, genotype, sep = "_"))
full_res_P635LR636Q_homvswt_log2FC_sig <- 
  full_res_P635LR636Q_homvswt_log2FC_sig %>% 
  mutate(sample_mutant_gen = paste(mutant, genotype, run, sep = "_"))
# Average read counts
full_res_P635LR636Q_homvswt_log2FC_sig <- 
  full_res_P635LR636Q_homvswt_log2FC_sig %>% 
  group_by(gene, mutant_gen) %>% 
  mutate(avg_value = mean(value))

# Average gene expression and log2FC
png("results/dotplot_STARalignment_P635LR636Q_hom_het_avg.png", 
    width = 8, height = 3.5, units = "in", res = 1000)
ggplot(full_res_P635LR636Q_homvswt_log2FC_sig, 
       aes(x = gene_name, y = mutant_gen)) +
  geom_point(aes(size = log1p(avg_value), fill = avg_log2FC), pch = 21) +
  scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white",
                       high = "red", space ="Lab" ) +
  labs(fill = "log2 Fold Change", size = "log gene expression") +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
dev.off()
# all samples
png("results/dotplot_STARalignment_P635LR636Q_hom_het_allsamples.png", 
    width = 8, height = 8, units = "in", res = 1000)
ggplot(full_res_P635LR636Q_homvswt_log2FC_sig, 
       aes(x = gene_name, y = sample_mutant_gen)) +
  geom_point(aes(size = log1p(value), fill = log2_FC), pch = 21) +
  scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white",
                       high = "red", space ="Lab" ) +
  labs(fill = "log2 Fold Change", size = "log gene expression") +
  xlab("gene") +
  ylab("sample") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
dev.off()

########### After editing ############

# Sample table
# Sample table
sample_table_after_base_editing <- read.table("/g/steinmetz/project/dcm_lgreads/mouse_bulkRNA/config/samples_afterbaseediting_metadata.txt",
                           header = TRUE)
# read.table(snakemake@input$meta, header = TRUE)
# Simplify sample table
sample_table_after_base_editing$sample <- paste("results/STAR", 
                             sample_table_after_base_editing$mutant, 
                             sample_table_after_base_editing$run, 
                             "Aligned.sortedByCoord.out.bam", sep = "/")

# Load data
featurecounts_P635L_after_base_editing <- readRDS("results/robject_featurecounts_P635L_after_base_editing.rds")
featurecounts_R636Q_after_base_editing <- readRDS("results/robject_featurecounts_R636Q_after_base_editing.rds")

# Combine counts
# Check if rownames
all(rownames(featurecounts_P635L_after_base_editing$counts) ==
  rownames(featurecounts_R636Q_after_base_editing$counts))
# Merge
featurecounts_after_base_editing <- 
  cbind(featurecounts_P635L_after_base_editing$counts,
       featurecounts_R636Q_after_base_editing$counts)

# PCA
# Prepare data for PCA
counts_after_base_editing <- t(featurecounts_after_base_editing)
counts_after_base_editing_wnames <- rownames_to_column(as.data.frame(counts_after_base_editing), 
                                             var = "sample")
counts_after_base_editing_wmeta <- left_join(counts_after_base_editing_wnames, 
                                             sample_table_after_base_editing, 
                                             keep = FALSE)
# Remove liver samples
counts_after_base_editing_wmeta_filtered <- dplyr::filter(counts_after_base_editing_wmeta, 
                                                          !is.na(treatment))
# Drop NA level
counts_after_base_editing_wmeta_filtered <- droplevels(counts_after_base_editing_wmeta_filtered)
# Recreate counts matrix
counts_after_base_editing_filtered <- dplyr::select(counts_after_base_editing_wmeta_filtered, contains("ENSMUSG"))

# Perform PCA
# All mice
pca_after_base_editing <- prcomp(counts_after_base_editing_filtered)
plot_pca_after_base_editing <- autoplot(pca_after_base_editing,
                              data = counts_after_base_editing_wmeta_filtered,
                              colour = "treatment",
                              shape = "mutant", 
                              size = 3,
                              alpha = 0.8, 
                              label = TRUE,  
                              label.label = "run",
                              label.position = position_nudge(x = 0, y = 0.015),
                              label.size = 2) +
  theme_bw() +
  ggtitle("PCA with all mice and all genes - after base editing experiment")
ggsave("results/plot_STARalignment_pca_after_base_editing.pdf",
       plot = plot_pca_after_base_editing,
       width = 7, height = 5, units = "in")

# DESeq2
# DESeqDataSetFromMatrix
# P635L
# Filter
sample_table_P635L_after_base_editing <- dplyr::filter(sample_table_after_base_editing, 
                                                       mutant == "P635L_after_base_editing")
# Extract counts
featurecounts_P635L_after_base_editing <- 
  featurecounts_P635L_after_base_editing$counts[,sample_table_P635L_after_base_editing$sample]
# Check order
if(all(colnames(featurecounts_P635L_after_base_editing) == sample_table_P635L_after_base_editing$sample)){
  colnames(featurecounts_P635L_after_base_editing) <- sample_table_P635L_after_base_editing$run
} else {
  # Rename columns
  stop(paste0("Check that the order of the columns of the counts for P635L mutant mice corresponds to the sample table"))
}
# Create deseq2 object
dds_P635L_after_base_editing <- DESeqDataSetFromMatrix(countData = featurecounts_P635L_after_base_editing,
                                    colData = sample_table_P635L_after_base_editing[,c("run", "treatment")],
                                    design = ~ treatment)
# Perform DE test
dds_P635L_after_base_editing <- DESeq(dds_P635L_after_base_editing)
# ALL comparisons
# P635L Nterm_SpRY_and_Cterm_SpRY_gRNA5 vs Nterm_NRTH_Abe8e_and_Cterm_gRNA5
res_P635L_after_base_editing_SpRYvsNRTH <- results(dds_P635L_after_base_editing, 
                                                contrast = c("treatment","Nterm_SpRY_and_Cterm_SpRY_gRNA5","Nterm_NRTH_Abe8e_and_Cterm_gRNA5"))
res_P635L_after_base_editing_SpRYvsNRTH_sig <- dplyr::filter(as.data.frame(res_P635L_after_base_editing_SpRYvsNRTH),
                                       padj < 0.05)
res_P635L_after_base_editing_SpRYvsNRTH_sig <- rownames_to_column(res_P635L_after_base_editing_SpRYvsNRTH_sig, var = "gene")
write.table(x = res_P635L_after_base_editing_SpRYvsNRTH_sig,
            file = "results/table_STARalignment_res_P635L_after_base_editing_SpRYvsNRTH_sig.txt",
            row.names = FALSE,
            col.names = TRUE)
res_P635L_after_base_editing_SpRYvsNRTH_sig_wgenename <- left_join(res_P635L_after_base_editing_SpRYvsNRTH_sig, 
                                                                  gen_vM29_annot_gene_uniq)
write.table(x = res_P635L_after_base_editing_SpRYvsNRTH_sig_wgenename,
            file = "results/table_STARalignment_res_P635L_after_base_editing_SpRYvsNRTH_sig_wgenename.txt",
            row.names = FALSE,
            col.names = TRUE)
# P635L Nterm_SpRY_and_Cterm_SpRY_gRNA5 vs PBS
res_P635L_after_base_editing_SpRYvsPBS <- results(dds_P635L_after_base_editing, 
                                                   contrast = c("treatment","Nterm_SpRY_and_Cterm_SpRY_gRNA5","PBS"))
res_P635L_after_base_editing_SpRYvsPBS_sig <- dplyr::filter(as.data.frame(res_P635L_after_base_editing_SpRYvsPBS),
                                                             padj < 0.05)
res_P635L_after_base_editing_SpRYvsPBS_sig <- rownames_to_column(res_P635L_after_base_editing_SpRYvsPBS_sig, var = "gene")
write.table(x = res_P635L_after_base_editing_SpRYvsPBS_sig,
            file = "results/table_STARalignment_res_P635L_after_base_editing_SpRYvsPBS_sig.txt",
            row.names = FALSE,
            col.names = TRUE)
res_P635L_after_base_editing_SpRYvsPBS_sig_wgenename <- left_join(res_P635L_after_base_editing_SpRYvsPBS_sig, 
                                                                   gen_vM29_annot_gene_uniq)
write.table(x = res_P635L_after_base_editing_SpRYvsPBS_sig_wgenename,
            file = "results/table_STARalignment_res_P635L_after_base_editing_SpRYvsPBS_sig_wgenename.txt",
            row.names = FALSE,
            col.names = TRUE)
# P635L Nterm_SpRY_and_Cterm_SpRY_gRNA5 vs none
res_P635L_after_base_editing_SpRYvsnone <- results(dds_P635L_after_base_editing, 
                                                  contrast = c("treatment","Nterm_SpRY_and_Cterm_SpRY_gRNA5","none"))
res_P635L_after_base_editing_SpRYvsnone_sig <- dplyr::filter(as.data.frame(res_P635L_after_base_editing_SpRYvsnone),
                                                            padj < 0.05)
res_P635L_after_base_editing_SpRYvsnone_sig <- rownames_to_column(res_P635L_after_base_editing_SpRYvsnone_sig, var = "gene")
write.table(x = res_P635L_after_base_editing_SpRYvsnone_sig,
            file = "results/table_STARalignment_res_P635L_after_base_editing_SpRYvsnone_sig.txt",
            row.names = FALSE,
            col.names = TRUE)
res_P635L_after_base_editing_SpRYvsnone_sig_wgenename <- left_join(res_P635L_after_base_editing_SpRYvsnone_sig, 
                                                                  gen_vM29_annot_gene_uniq)
write.table(x = res_P635L_after_base_editing_SpRYvsnone_sig_wgenename,
            file = "results/table_STARalignment_res_P635L_after_base_editing_SpRYvsnone_sig_wgenename.txt",
            row.names = FALSE,
            col.names = TRUE)
# P635L Nterm_NRTH_Abe8e_and_Cterm_gRNA5 vs PBS
res_P635L_after_base_editing_NRTHvsPBS <- results(dds_P635L_after_base_editing, 
                                                  contrast = c("treatment","Nterm_NRTH_Abe8e_and_Cterm_gRNA5","PBS"))
res_P635L_after_base_editing_NRTHvsPBS_sig <- dplyr::filter(as.data.frame(res_P635L_after_base_editing_NRTHvsPBS),
                                                            padj < 0.05)
res_P635L_after_base_editing_NRTHvsPBS_sig <- rownames_to_column(res_P635L_after_base_editing_NRTHvsPBS_sig, var = "gene")
write.table(x = res_P635L_after_base_editing_NRTHvsPBS_sig,
            file = "results/table_STARalignment_res_P635L_after_base_editing_NRTHvsPBS_sig.txt",
            row.names = FALSE,
            col.names = TRUE)
res_P635L_after_base_editing_NRTHvsPBS_sig_wgenename <- left_join(res_P635L_after_base_editing_NRTHvsPBS_sig, 
                                                                   gen_vM29_annot_gene_uniq)
write.table(x = res_P635L_after_base_editing_NRTHvsPBS_sig_wgenename,
            file = "results/table_STARalignment_res_P635L_after_base_editing_NRTHvsPBS_sig_wgenename.txt",
            row.names = FALSE,
            col.names = TRUE)
# P635L Nterm_NRTH_Abe8e_and_Cterm_gRNA5 vs none
res_P635L_after_base_editing_NRTHvsnone <- results(dds_P635L_after_base_editing, 
                                                  contrast = c("treatment","Nterm_NRTH_Abe8e_and_Cterm_gRNA5","none"))
res_P635L_after_base_editing_NRTHvsnone_sig <- dplyr::filter(as.data.frame(res_P635L_after_base_editing_NRTHvsnone),
                                                            padj < 0.05)
res_P635L_after_base_editing_NRTHvsnone_sig <- rownames_to_column(res_P635L_after_base_editing_NRTHvsnone_sig, var = "gene")
write.table(x = res_P635L_after_base_editing_NRTHvsnone_sig,
            file = "results/table_STARalignment_res_P635L_after_base_editing_NRTHvsnone_sig.txt",
            row.names = FALSE,
            col.names = TRUE)
res_P635L_after_base_editing_NRTHvsnone_sig_wgenename <- left_join(res_P635L_after_base_editing_NRTHvsnone_sig, 
                                                                  gen_vM29_annot_gene_uniq)
write.table(x = res_P635L_after_base_editing_NRTHvsnone_sig_wgenename,
            file = "results/table_STARalignment_res_P635L_after_base_editing_NRTHvsnone_sig_wgenename.txt",
            row.names = FALSE,
            col.names = TRUE)
# P635L PBS vs none
res_P635L_after_base_editing_PBSvsnone <- results(dds_P635L_after_base_editing, 
                                                   contrast = c("treatment","PBS","none"))
res_P635L_after_base_editing_PBSvsnone_sig <- dplyr::filter(as.data.frame(res_P635L_after_base_editing_PBSvsnone),
                                                             padj < 0.05)
res_P635L_after_base_editing_PBSvsnone_sig <- rownames_to_column(res_P635L_after_base_editing_PBSvsnone_sig, var = "gene")
write.table(x = res_P635L_after_base_editing_PBSvsnone_sig,
            file = "results/table_STARalignment_res_P635L_after_base_editing_PBSvsnone_sig.txt",
            row.names = FALSE,
            col.names = TRUE)
res_P635L_after_base_editing_PBSvsnone_sig_wgenename <- left_join(res_P635L_after_base_editing_PBSvsnone_sig, 
                                                                  gen_vM29_annot_gene_uniq)
write.table(x = res_P635L_after_base_editing_PBSvsnone_sig_wgenename,
            file = "results/table_STARalignment_res_P635L_after_base_editing_PBSvsnone_sig_wgenename.txt",
            row.names = FALSE,
            col.names = TRUE)

# R636Q
# Filter
sample_table_R636Q_after_base_editing <- dplyr::filter(sample_table_after_base_editing, 
                                                       mutant == "R636Q_after_base_editing")
# Extract counts
featurecounts_R636Q_after_base_editing <- 
  featurecounts_R636Q_after_base_editing$counts[,sample_table_R636Q_after_base_editing$sample]
# Check order
if(all(colnames(featurecounts_R636Q_after_base_editing) == sample_table_R636Q_after_base_editing$sample)){
  colnames(featurecounts_R636Q_after_base_editing) <- sample_table_R636Q_after_base_editing$run
} else {
  # Rename columns
  stop(paste0("Check that the order of the columns of the counts for R636Q mutant mice corresponds to the sample table"))
}
# Create deseq2 object
dds_R636Q_after_base_editing <- DESeqDataSetFromMatrix(countData = featurecounts_R636Q_after_base_editing,
                                                       colData = sample_table_R636Q_after_base_editing[,c("run", "treatment")],
                                                       design = ~ treatment)
# Perform DE test
dds_R636Q_after_base_editing <- DESeq(dds_R636Q_after_base_editing)
# ALL comparisons
# R636Q Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6 vs PBS
res_R636Q_after_base_editing_NRTHvsPBS <- results(dds_R636Q_after_base_editing, 
                                                   contrast = c("treatment","Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6","PBS"))
res_R636Q_after_base_editing_NRTHvsPBS_sig <- dplyr::filter(as.data.frame(res_R636Q_after_base_editing_NRTHvsPBS),
                                                             padj < 0.05)
res_R636Q_after_base_editing_NRTHvsPBS_sig <- rownames_to_column(res_R636Q_after_base_editing_NRTHvsPBS_sig, var = "gene")
write.table(x = res_R636Q_after_base_editing_NRTHvsPBS_sig,
            file = "results/table_STARalignment_res_R636Q_after_base_editing_NRTHvsPBS_sig.txt",
            row.names = FALSE,
            col.names = TRUE)
res_R636Q_after_base_editing_NRTHvsPBS_sig_wgenename <- left_join(res_R636Q_after_base_editing_NRTHvsPBS_sig, 
                                                                   gen_vM29_annot_gene_uniq)
write.table(x = res_R636Q_after_base_editing_NRTHvsPBS_sig_wgenename,
            file = "results/table_STARalignment_res_R636Q_after_base_editing_NRTHvsPBS_sig_wgenename.txt",
            row.names = FALSE,
            col.names = TRUE)
# R636Q Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6 vs none
res_R636Q_after_base_editing_NRTHvsnone <- results(dds_R636Q_after_base_editing, 
                                                  contrast = c("treatment","Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6","none"))
res_R636Q_after_base_editing_NRTHvsnone_sig <- dplyr::filter(as.data.frame(res_R636Q_after_base_editing_NRTHvsnone),
                                                            padj < 0.05)
res_R636Q_after_base_editing_NRTHvsnone_sig <- rownames_to_column(res_R636Q_after_base_editing_NRTHvsnone_sig, var = "gene")
write.table(x = res_R636Q_after_base_editing_NRTHvsnone_sig,
            file = "results/table_STARalignment_res_R636Q_after_base_editing_NRTHvsnone_sig.txt",
            row.names = FALSE,
            col.names = TRUE)
res_R636Q_after_base_editing_NRTHvsnone_sig_wgenename <- left_join(res_R636Q_after_base_editing_NRTHvsnone_sig, 
                                                                  gen_vM29_annot_gene_uniq)
write.table(x = res_R636Q_after_base_editing_NRTHvsnone_sig_wgenename,
            file = "results/table_STARalignment_res_R636Q_after_base_editing_NRTHvsnone_sig_wgenename.txt",
            row.names = FALSE,
            col.names = TRUE)
# R636Q PBS vs none
res_R636Q_after_base_editing_PBSvsnone <- results(dds_R636Q_after_base_editing, 
                                                   contrast = c("treatment","PBS","none"))
res_R636Q_after_base_editing_PBSvsnone_sig <- dplyr::filter(as.data.frame(res_R636Q_after_base_editing_PBSvsnone),
                                                             padj < 0.05)
res_R636Q_after_base_editing_PBSvsnone_sig <- rownames_to_column(res_R636Q_after_base_editing_PBSvsnone_sig, 
                                                                 var = "gene")
write.table(x = res_R636Q_after_base_editing_PBSvsnone_sig,
            file = "results/table_STARalignment_res_R636Q_after_base_editing_PBSvsnone_sig.txt",
            row.names = FALSE,
            col.names = TRUE)
res_R636Q_after_base_editing_PBSvsnone_sig_wgenename <- left_join(res_R636Q_after_base_editing_PBSvsnone_sig, 
          gen_vM29_annot_gene_uniq)
write.table(x = res_R636Q_after_base_editing_PBSvsnone_sig_wgenename,
            file = "results/table_STARalignment_res_R636Q_after_base_editing_PBSvsnone_sig_wgenename.txt",
            row.names = FALSE,
            col.names = TRUE)

# Heatmaps
# P635L
pdf("results/heatmap_STARalignment_P635L_after_base_editing_SpRYvsNRTH_top20.pdf")
function_my_heatmap(mutation = "P635L",
                    dds = dds_P635L_after_base_editing,
                    res = res_P635L_after_base_editing_SpRYvsNRTH_sig,
                    sample_tab = sample_table_P635L_after_base_editing, 
                    annot = gen_vM29_annot_gene_uniq,
                    gen = c("Nterm_SpRY_and_Cterm_SpRY_gRNA5", "Nterm_NRTH_Abe8e_and_Cterm_gRNA5"),
                    filter_col = "treatment",
                    top = 20)
dev.off()
pdf("results/heatmap_STARalignment_P635L_after_base_editing_SpRYvsPBS_top20.pdf")
function_my_heatmap(mutation = "P635L",
                    dds = dds_P635L_after_base_editing,
                    res = res_P635L_after_base_editing_SpRYvsPBS_sig,
                    sample_tab = sample_table_P635L_after_base_editing, 
                    annot = gen_vM29_annot_gene_uniq,
                    gen = c("Nterm_SpRY_and_Cterm_SpRY_gRNA5", "PBS"),
                    filter_col = "treatment",
                    top = 20)
dev.off()
pdf("results/heatmap_STARalignment_P635L_after_base_editing_SpRYvsnone_top20.pdf")
function_my_heatmap(mutation = "P635L",
                    dds = dds_P635L_after_base_editing,
                    res = res_P635L_after_base_editing_SpRYvsnone_sig,
                    sample_tab = sample_table_P635L_after_base_editing, 
                    annot = gen_vM29_annot_gene_uniq,
                    gen = c("Nterm_SpRY_and_Cterm_SpRY_gRNA5", "none"),
                    filter_col = "treatment",
                    top = 20)
dev.off()
pdf("results/heatmap_STARalignment_P635L_after_base_editing_NRTHvsPBS_top20.pdf")
function_my_heatmap(mutation = "P635L",
                    dds = dds_P635L_after_base_editing,
                    res = res_P635L_after_base_editing_NRTHvsPBS_sig,
                    sample_tab = sample_table_P635L_after_base_editing, 
                    annot = gen_vM29_annot_gene_uniq,
                    gen = c("Nterm_NRTH_Abe8e_and_Cterm_gRNA5", "PBS"),
                    filter_col = "treatment",
                    top = 20)
dev.off()
pdf("results/heatmap_STARalignment_P635L_after_base_editing_NRTHvsnone_top20.pdf")
function_my_heatmap(mutation = "P635L",
                    dds = dds_P635L_after_base_editing,
                    res = res_P635L_after_base_editing_NRTHvsnone_sig,
                    sample_tab = sample_table_P635L_after_base_editing, 
                    annot = gen_vM29_annot_gene_uniq,
                    gen = c("Nterm_NRTH_Abe8e_and_Cterm_gRNA5", "none"),
                    filter_col = "treatment",
                    top = 20)
dev.off()
pdf("results/heatmap_STARalignment_P635L_after_base_editing_PBSvsnone_top20.pdf")
function_my_heatmap(mutation = "P635L",
                    dds = dds_P635L_after_base_editing,
                    res = res_P635L_after_base_editing_PBSvsnone_sig,
                    sample_tab = sample_table_P635L_after_base_editing, 
                    annot = gen_vM29_annot_gene_uniq,
                    gen = c("PBS", "none"),
                    filter_col = "treatment",
                    top = 20)
dev.off()
# R636Q
pdf("results/heatmap_STARalignment_R636Q_after_base_editing_PBSvsnone_top20.pdf")
function_my_heatmap(mutation = "R636Q",
                    dds = dds_R636Q_after_base_editing,
                    res = res_R636Q_after_base_editing_PBSvsnone_sig,
                    sample_tab = sample_table_R636Q_after_base_editing, 
                    annot = gen_vM29_annot_gene_uniq,
                    gen = c("PBS", "none"),
                    filter_col = "treatment",
                    top = 20)
dev.off()
pdf("results/heatmap_STARalignment_R636Q_after_base_editing_NRTHvsnone_top20.pdf")
function_my_heatmap(mutation = "R636Q",
                    dds = dds_R636Q_after_base_editing,
                    res = res_R636Q_after_base_editing_NRTHvsnone_sig,
                    sample_tab = sample_table_R636Q_after_base_editing, 
                    annot = gen_vM29_annot_gene_uniq,
                    gen = c("Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6", "none"),
                    filter_col = "treatment",
                    top = 20)
dev.off()
pdf("results/heatmap_STARalignment_R636Q_after_base_editing_NRTHvsPBS_top20.pdf")
function_my_heatmap(mutation = "R636Q",
                    dds = dds_R636Q_after_base_editing,
                    res = res_R636Q_after_base_editing_NRTHvsPBS_sig,
                    sample_tab = sample_table_R636Q_after_base_editing, 
                    annot = gen_vM29_annot_gene_uniq,
                    gen = c("Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6", "PBS"),
                    filter_col = "treatment",
                    top = 20)
dev.off()

# Merged heatmap
# Assays
merged_assay_geneexp_after_base_editing <- cbind(counts(dds_P635L_after_base_editing, normalized=TRUE), 
                                                 counts(dds_R636Q_after_base_editing, normalized=TRUE))
# Sample tables
merged_sample_table_after_base_editing <- rbind(sample_table_P635L_after_base_editing,
                             sample_table_R636Q_after_base_editing)
simpl_merged_sample_table_after_base_editing <- dplyr::select(merged_sample_table_after_base_editing,
                                                              run, mutant, treatment)
row.names(simpl_merged_sample_table_after_base_editing) <- 
  simpl_merged_sample_table_after_base_editing$run
# Plot
pdf("results/heatmap_STARalignment_P635L_and_R636Q_after_base_editing_allgenes.pdf",
    width = 10)
pheatmap(merged_assay_geneexp_after_base_editing,
         cluster_rows = FALSE, 
         show_rownames = FALSE,
         annotation_col = simpl_merged_sample_table_after_base_editing[,c("mutant", "treatment")],
         scale = "row")
dev.off()

# List of gene provided by Markus 09.13.2022
# Pulled out from the initial experiments
# Common genes
DEgenelist_initialexp_common <- read.table("results/DEgenelist_initialexp_common_p_inf_1e-10.txt",
                                           header = TRUE)
# Filter
merged_assay_geneexp_after_base_editing_tofiltered <- rownames_to_column(as.data.frame(merged_assay_geneexp_after_base_editing), var = "gene")
assay_after_base_editing_DEgenelist_initialexp_common <- dplyr::left_join(DEgenelist_initialexp_common, 
                 merged_assay_geneexp_after_base_editing_tofiltered)
assay_after_base_editing_DEgenelist_initialexp_common <- column_to_rownames(assay_after_base_editing_DEgenelist_initialexp_common,
                   var = "gene_name")
assay_after_base_editing_DEgenelist_initialexp_common <- dplyr::select(assay_after_base_editing_DEgenelist_initialexp_common, !(gene))
# Plot
pdf("results/heatmap_STARalignment_P635L_and_R636Q_after_base_editing_DEgenelist_initialexp_common_p_inf_1e-10.pdf",
    width = 10, height = 20)
pheatmap(assay_after_base_editing_DEgenelist_initialexp_common,
         cluster_rows = TRUE, 
         show_rownames = TRUE,
         annotation_col = simpl_merged_sample_table_after_base_editing[,c("mutant", "treatment")],
         scale = "row",
         fontsize = 5,
         main = "P635L and R636Q after base editing - Common DE genes from initial experiment p < 1e-10")
dev.off()
# P635L genes
DEgenelist_initialexp_P635L <- read.table("results/DEgenelist_initialexp_P635L_p_inf_1e-10.txt",
                                           header = TRUE)
# Filter
P635L_merged_assay_geneexp_after_base_editing_tofiltered <- rownames_to_column(as.data.frame(counts(dds_P635L_after_base_editing, normalized=TRUE)), var = "gene")
assay_after_base_editing_DEgenelist_initialexp_P635L <- dplyr::left_join(DEgenelist_initialexp_P635L, 
                                                                         P635L_merged_assay_geneexp_after_base_editing_tofiltered)
assay_after_base_editing_DEgenelist_initialexp_P635L <- column_to_rownames(assay_after_base_editing_DEgenelist_initialexp_P635L,
                                                                            var = "gene_name")
assay_after_base_editing_DEgenelist_initialexp_P635L <- dplyr::select(assay_after_base_editing_DEgenelist_initialexp_P635L, !(gene))
# Plot
pdf("results/heatmap_STARalignment_P635L_after_base_editing_DEgenelist_initialexp_P635L_p_inf_1e-10.pdf",
    width = 10, height = 20)
pheatmap(assay_after_base_editing_DEgenelist_initialexp_P635L,
         cluster_rows = TRUE, 
         show_rownames = TRUE,
         annotation_col = simpl_merged_sample_table_after_base_editing[,c("mutant", "treatment")],
         scale = "row",
         fontsize = 5,
         main = "P635L after base editing - DE genes from P635L initial experiment p < 1e-10")
dev.off()
# R636Q genes
DEgenelist_initialexp_R636Q <- read.table("results/DEgenelist_initialexp_R636Q_p_inf_1e-10.txt",
                                          header = TRUE)
# Filter
R636Q_merged_assay_geneexp_after_base_editing_tofiltered <- rownames_to_column(as.data.frame(counts(dds_R636Q_after_base_editing, normalized=TRUE)), var = "gene")
assay_after_base_editing_DEgenelist_initialexp_R636Q <- dplyr::left_join(DEgenelist_initialexp_R636Q, 
                                                                         R636Q_merged_assay_geneexp_after_base_editing_tofiltered)
assay_after_base_editing_DEgenelist_initialexp_R636Q <- column_to_rownames(assay_after_base_editing_DEgenelist_initialexp_R636Q,
                                                                           var = "gene_name")
assay_after_base_editing_DEgenelist_initialexp_R636Q <- dplyr::select(assay_after_base_editing_DEgenelist_initialexp_R636Q, !(gene))
# Plot
pdf("results/heatmap_STARalignment_R636Q_after_base_editing_DEgenelist_initialexp_R636Q_p_inf_1e-10.pdf",
    width = 10, height = 20)
pheatmap(assay_after_base_editing_DEgenelist_initialexp_R636Q,
         cluster_rows = TRUE, 
         show_rownames = TRUE,
         annotation_col = simpl_merged_sample_table_after_base_editing[,c("mutant", "treatment")],
         scale = "row",
         fontsize = 5,
         main = "R636Q after base editing - DE genes from R636Q initial experiment p < 1e-10")
dev.off()
