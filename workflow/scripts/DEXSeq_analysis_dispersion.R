# DEXSeq estimate dispersion

# Load libraries
library(DEXSeq)
library(tidyverse)
library(BiocParallel)

# Mutation
mut <- snakemake@wildcards$mutation

# Read counts file from STAR alignment
# List files
count_files <- list.files(paste0("./results/STAR/", mut), 
                          pattern="read_counts_sorted_reverse_clean.txt$", 
                          recursive = TRUE,
                          full.names = TRUE)
# Annotation file
flattened_file <- list.files("./results/gencode_GRCm39",
                             pattern = "*.gff", 
                             full.names = TRUE)
# Sample table
sample_table <- read.table(snakemake@input$meta, header = TRUE)
# sample_table <- read.table("/g/steinmetz/project/dcm_lgreads/mouse_bulkRNA/config/samples_metadata.txt",
#                            header = TRUE)
# Filter
sample_table <- dplyr::filter(sample_table, mutant == mut)
# Add library information
sample_table$libType <- "paired-end" 
# Simplify sample table
sample_table$rowname_IDs <- paste(sample_table$mutant, sample_table$run, "read_counts_sorted_reverse_clean", sep = "/")
sample_table <- column_to_rownames(sample_table, var = "rowname_IDs")
sample_table_dex <- sample_table[,c("genotype", "libType")]
sample_table_dex <- dplyr::rename(sample_table_dex, condition = genotype)

# 3 conditions
# Create DEXSeq object
dxd <- DEXSeqDataSetFromHTSeq(
  count_files,
  sampleData = sample_table_dex,
  design= ~ sample + exon + condition:exon,
  flattenedfile = flattened_file)
# the first 15 (we have 15 samples) corresponding to the number of reads mapping 
# to out exonic regions and the last 15 correspond to the sum of the counts
# mapping to the rest of the exons from the same gene on each sample.
split(seq_len(ncol(dxd)), colData(dxd)$exon)
# rows are labelled with gene IDs , followed by a colon and the counting bin number
head(DEXSeq::featureCounts(dxd), 5)
# details on the counting bins
head(rowRanges(dxd), 3)
# Design table
sampleAnnotation(dxd)

# Normalization
# DEXSeq uses the same method as DESeq2, which is provided in the function estimateSizeFactors().
dxd <- estimateSizeFactors(dxd)

# Dispersion estimation
# To test for differential exon usage, we need to estimate the variability of the data.
# The goal of this step is to distinguish technical and biological variation (noise) 
# from real effects on exon usage due to the different conditions
# TO NOTE: samples with the same experimental condition are considered as replicates 
dxd <- estimateDispersions(dxd, BPPARAM = MulticoreParam(snakemake@threads))
# Save object
save(dxd, file = paste0("results/DEXSeq/robject_", {mut}, "_dexseq_estdisp_3conditions.rda"))
# plot
png(paste0("results/DEXSeq/plot_", mut, "_dxd_dispersion_estimate_3conditions.png"))
plotDispEsts(dxd)
dev.off()

# hom vs WT
# Filter
sample_table_dex_wt_hom <- dplyr::filter(sample_table_dex, condition != "het")
idx_count_files_wt_hom <- count_files %in% paste0("./results/STAR/", rownames(sample_table_dex_wt_hom), ".txt")
count_files_wt_hom <- count_files[idx_count_files_wt_hom]
# Create DEXSeq object
dxd_wt_hom <- DEXSeqDataSetFromHTSeq(
  count_files_wt_hom,
  sampleData = sample_table_dex_wt_hom,
  design= ~ sample + exon + condition:exon,
  flattenedfile = flattened_file)
split(seq_len(ncol(dxd_wt_hom)), colData(dxd_wt_hom)$exon)
# Normalization
dxd_wt_hom <- estimateSizeFactors(dxd_wt_hom)
# Dispersion estimation
dxd_wt_hom <- estimateDispersions(dxd_wt_hom, BPPARAM = MulticoreParam(snakemake@threads))
# Save object
save(dxd_wt_hom, file = paste0("results/DEXSeq/robject_", {mut}, "_dexseq_estdisp_WTvsHOM.rda"))
# plot
png(paste0("results/DEXSeq/plot_", mut, "_dxd_dispersion_estimate_WTvsHOM.png"))
plotDispEsts(dxd_wt_hom)
dev.off()
