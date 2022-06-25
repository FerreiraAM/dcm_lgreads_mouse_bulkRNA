# DEBUG SCRIPT
# save.image("deseq_analysis.rda")
# stop()

set.seed(1234)

# Load libraries:
library(tidyverse) 
library(DESeq2)
library(tximport)
library(rtracklayer)
library(ggplot2)
library(ggfortify)
library(ggrepel)
library(pheatmap)

# Data analysis 

# Gene ID data frame
tx_g_ids <- read.table("resources/mus_musculus/transcripts_to_genes.txt")
tx_g_ids <- tx_g_ids[,c(1,3)]
colnames(tx_g_ids) <- c("TXNAME", "GENEID")
# unique(tx_g_ids)
# Metadata
# Read samples metadata file
samples_metadata <- read.table(snakemake@input$meta, header = TRUE)

# Import data from Kallisto alignment
# Files list
files <- file.path(samples_metadata$results_dir, 
                   samples_metadata$mutant, 
                   samples_metadata$run, 
                   "abundance.tsv")
names(files) <- samples_metadata$run
# Import
all_txi <- tximport(files, 
                    type = "kallisto", 
                    tx2gene = unique(tx_g_ids))
# ----------------------------------------------------------------------------
# # 2829 missing transcripts from tx2gene
# checktx <- read.table(files[1], sep = "\t", header = TRUE)
# intx2gene <- checktx$target_id %in% tx_g_ids$TXNAME
# notintx2gene <- which(!intx2gene)
# txnotintx2gene <- checktx$target_id[notintx2gene]
# length(txnotintx2gene) # 2829
# # Remove version
# first_word <- function(mystring){
#           unlist(strsplit(mystring, "[.]"))[1]
# }
# txnotintx2genenotv <- sapply(txnotintx2gene, first_word)
# # Download gtf file
# gtf_file <- import("/g/steinmetz/project/dcm_lgreads/mouse_bulkRNA/resources/mus_musculus/Mus_musculus.GRCm38.96.gtf", 
#        format = "gtf")
# all(txnotintx2genenotv %in% mcols(gtf_file)$transcript_id) #not all of them are in the gtf file
# gtf_txnotintx2genenotv <- dplyr::filter(as.data.frame(mcols(gtf_file)), 
#                                         transcript_id %in% txnotintx2genenotv)
# unique(gtf_txnotintx2genenotv$gene_name)
# ----------------------------------------------------------------------------
# Number of reads
# head(all_txi$counts)
total_nb_reads <- colSums(all_txi$counts)
total_nb_reads <- rownames_to_column(as.data.frame(total_nb_reads))
colnames(total_nb_reads)[1] <- c("run")
total_nb_reads_merge_w_metadata <- left_join(samples_metadata, total_nb_reads)
write.csv(total_nb_reads_merge_w_metadata[,c("mutant", "run", "genotype", "total_nb_reads")], 
          file = "results/total_nb_reads.csv", 
          col.names = TRUE, row.names = FALSE)

# PCA
# Prepare data for PCA
counts_all_mice <- t(all_txi$counts)
counts_all_mice_wnames <- rownames_to_column(as.data.frame(counts_all_mice), var = "run")
counts_all_mice_wmeta <- left_join(counts_all_mice_wnames, samples_metadata)
# Perform PCA
# All mice
pca_all_mice <- prcomp(counts_all_mice)
plot_pca_all_mice <- autoplot(pca_all_mice, data = counts_all_mice_wmeta, colour = "mutant") +
  theme_bw() +
  ggtitle("PCA with all mice and all genes")
ggsave("results/plot_pca_all_mice_genes.pdf", 
       plot = plot_pca_all_mice)
# P635L
# Filter
counts_P635L_mice_wmeta <- dplyr::filter(counts_all_mice_wmeta, 
                                         mutant == "P635L")
counts_P635L_mice <- dplyr::select(counts_P635L_mice_wmeta, 
                                   !c("run", "mutant", "genotype", "results_dir"))
pca_P635L_mice <- prcomp(counts_P635L_mice)
plot_pca_P635L_mice <- autoplot(pca_P635L_mice, 
                     data = counts_P635L_mice_wmeta, 
                     colour = "genotype") +
  geom_text_repel(aes(label = run, color = genotype)) +
  theme_bw() +
  ggtitle("PCA with P635L mice and all genes")
ggsave("results/plot_pca_P635L_mice_genes.pdf", 
       plot = plot_pca_P635L_mice)
# R636Q
# Filter
counts_R636Q_mice_wmeta <- dplyr::filter(counts_all_mice_wmeta, 
                                         mutant == "R636Q")
counts_R636Q_mice <- dplyr::select(counts_R636Q_mice_wmeta, 
                                   !c("run", "mutant", "genotype", "results_dir"))
pca_R636Q_mice <- prcomp(counts_R636Q_mice)
plot_pca_R636Q_mice <- autoplot(pca_R636Q_mice, 
                                data = counts_R636Q_mice_wmeta, 
                                colour = "genotype") +
  geom_text_repel(aes(label = run, color = genotype)) +
  theme_bw() +
  ggtitle("PCA with R636Q mice and all genes")
ggsave("results/plot_pca_R636Q_mice_genes.pdf", 
       plot = plot_pca_R636Q_mice)

# DESeq2
# P635L 
# Filter Metdata
samples_metadata_P635L <- dplyr::filter(samples_metadata, mutant == "P635L")
sample_table_P635L <- data.frame(genotype = samples_metadata_P635L$genotype,
                                 row.names = samples_metadata_P635L$run)
# Filter transcript list
fct_sel <- function(data, vect){
  df <- as.data.frame(data)
  sub <- dplyr::select(df, vect)
  as.matrix(sub)
}
ls_txi_P635L <- lapply(all_txi[c("abundance", "counts", "length")], 
       FUN = fct_sel, 
       vect = samples_metadata_P635L$run)
ls_txi_P635L$countsFromAbundance <- all_txi$countsFromAbundance
# Create DESeq object
dds_P635L <- DESeqDataSetFromTximport(txi = ls_txi_P635L, sample_table_P635L, ~genotype)
# Differential test
dds_P635L <- DESeq(dds_P635L)
# P635L hom vs WT
res_P635L_homvswt <- results(dds_P635L, contrast=c("genotype","hom","wt"))
res_P635L_homvswt_sig <- dplyr::filter(as.data.frame(res_P635L_homvswt), 
                                       padj < 0.05)
res_P635L_homvswt_sig <- rownames_to_column(res_P635L_homvswt_sig, var = "run")
write.table(x = res_P635L_homvswt_sig,
            file = "results/table_res_P635L_homvswt_sig.txt",
            row.names = FALSE,
            col.names = TRUE)
# P635L het vs WT
res_P635L_hetvswt <- results(dds_P635L, contrast=c("genotype","het","wt"))
res_P635L_hetvswt_sig <- dplyr::filter(as.data.frame(res_P635L_hetvswt), 
                                       padj < 0.05)
res_P635L_hetvswt_sig <- rownames_to_column(res_P635L_hetvswt_sig, var = "run")
write.table(x = res_P635L_hetvswt_sig,
            file = "results/table_res_P635L_hetvswt_sig.txt",
            row.names = FALSE,
            col.names = TRUE)
# R636Q 
# Filter Metdata
samples_metadata_R636Q <- dplyr::filter(samples_metadata, mutant == "R636Q")
sample_table_R636Q <- data.frame(genotype = samples_metadata_R636Q$genotype,
                                 row.names = samples_metadata_R636Q$run)
# Filter transcript list
ls_txi_R636Q <- lapply(all_txi[c("abundance", "counts", "length")], 
                       FUN = fct_sel, 
                       vect = samples_metadata_R636Q$run)
ls_txi_R636Q$countsFromAbundance <- all_txi$countsFromAbundance
# Create DESeq object
dds_R636Q <- DESeqDataSetFromTximport(txi = ls_txi_R636Q, sample_table_R636Q, ~genotype)
# Differential test
dds_R636Q <- DESeq(dds_R636Q)
# R636Q hom vs WT
res_R636Q_homvswt <- results(dds_R636Q, contrast=c("genotype","hom","wt"))
res_R636Q_homvswt_sig <- dplyr::filter(as.data.frame(res_R636Q_homvswt), 
                                       padj < 0.05)
res_R636Q_homvswt_sig <- rownames_to_column(res_R636Q_homvswt_sig, var = "run")
write.table(x = res_R636Q_homvswt_sig,
            file = "results/table_res_R636Q_homvswt_sig.txt",
            row.names = FALSE,
            col.names = TRUE)
# R636Q het vs WT
res_R636Q_hetvswt <- results(dds_R636Q, contrast=c("genotype","het","wt"))
res_R636Q_hetvswt_sig <- dplyr::filter(as.data.frame(res_R636Q_hetvswt), 
                                       padj < 0.05)
res_R636Q_hetvswt_sig <- rownames_to_column(res_R636Q_hetvswt_sig, var = "run")
write.table(x = res_R636Q_hetvswt_sig,
            file = "results/table_res_R636Q_hetvswt_sig.txt",
            row.names = FALSE,
            col.names = TRUE)

# Expression heat map: R636Q hom, het, P635L hom, het, WT
# R636Q 
# Normalize the data
ntd_R636Q <- normTransform(dds_R636Q)
# Select significant genes
select_R636Q_hom <- unique(c(res_R636Q_hetvswt_sig$run, res_R636Q_homvswt_sig$run))
# select_R636Q_hom <- order(rowMeans(counts(dds_R636Q, normalized=TRUE)),
                # decreasing=TRUE)[1:20]
df_R636Q_hom <- as.data.frame(colData(dds_R636Q))
df_R636Q_hom$genotype <- factor(df_R636Q_hom$genotype, ordered = TRUE)
pdf(file = "results/heatmap_R636Q.pdf")
pheatmap(assay(ntd_R636Q)[select_R636Q_hom,], 
         cluster_rows=FALSE, 
         show_rownames=FALSE,
         cluster_cols=FALSE, 
         annotation_col=df_R636Q_hom)
dev.off()
# P635L 
# Normalize the data
ntd_P635L <- normTransform(dds_P635L)
# Select significant genes
select_P635L_hom <- unique(c(res_P635L_hetvswt_sig$run, res_P635L_homvswt_sig$run))
# select_P635L_hom <- order(rowMeans(counts(dds_P635L, normalized=TRUE)),
# decreasing=TRUE)[1:20]
df_P635L_hom <- as.data.frame(colData(dds_P635L))
df_P635L_hom$genotype <- factor(df_P635L_hom$genotype, ordered = TRUE)
pdf(file = "results/heatmap_P635L.pdf")
pheatmap(assay(ntd_P635L)[select_P635L_hom,], 
         cluster_rows=FALSE, 
         show_rownames=FALSE,
         cluster_cols=FALSE, 
         annotation_col=df_P635L_hom)
dev.off()
# select <- order(rowMeans(counts(dds,normalized=TRUE)),
#                 decreasing=TRUE)[1:20]
# df <- as.data.frame(colData(dds)[,c("condition","type")])
# pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
#          cluster_cols=FALSE, annotation_col=df)