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

# # Order genes based on adjusted p-values
# res_R636Q_homvswt_sig_ordered <- res_R636Q_homvswt_sig[order(res_R636Q_homvswt_sig$padj),]
# # Select hom and wt samples
# df_hom_wt_R636Q_table <- dplyr::filter(sample_table_R636Q, genotype == "hom" | genotype == "wt") %>% 
#   select(run, genotype)
# row.names(df_hom_wt_R636Q_table) <- df_hom_wt_R636Q_table$run
# # Extract top 20 genes
# counts_R636Q_homvswt_sig_allsamples <- counts(dds_R636Q, normalized=TRUE)[res_R636Q_homvswt_sig_ordered[1:20,"gene"],]
# counts_R636Q_homvswt_sig_allsamples <- rownames_to_column(as.data.frame(counts_R636Q_homvswt_sig_allsamples), 
#                                                var = "gene")
# # Select hom and wt samples
# counts_R636Q_homvswt_sig <- dplyr::select(counts_R636Q_homvswt_sig_allsamples, 
#                                           c(df_hom_wt_R636Q_table$run, "gene"))
# # Add gene name
# counts_R636Q_homvswt_sig_wgname <- left_join(counts_R636Q_homvswt_sig, 
#                                              gen_vM29_annot_gene_uniq)
# counts_R636Q_homvswt_sig_wgname <- column_to_rownames(counts_R636Q_homvswt_sig_wgname, "gene_name")
# # Plot
# pdf("results/heatmap_R636Q_homvswt_top20.pdf")
# pheatmap(counts_R636Q_homvswt_sig_wgname[,1:length(df_hom_wt_R636Q_table$run)], 
#          cluster_rows = TRUE, 
#          show_rownames = TRUE,
#          show_colnames = FALSE,
#          cluster_cols = TRUE,
#          annotation_col = df_hom_wt_R636Q_table,
#          scale = "row")
# dev.off()