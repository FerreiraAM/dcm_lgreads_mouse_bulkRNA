# DEXSeq analysis

# Load libraries
library(DEXSeq)
# library(tximport)
library(tidyverse)
library(rtracklayer)

# Read counts file from STAR alignment
# List files
count_files <- list.files("./results/STAR", 
                          pattern="read_counts_sorted_clean.txt$", 
                          recursive = TRUE,
                          full.names = TRUE)
# Annotation file
flattened_file <- list.files("./results/gencode_GRCm39",
                             pattern = "*.gff", 
                             full.names = TRUE)
# Sample table
sample_table <- read.table("/g/steinmetz/project/dcm_lgreads/mouse_bulkRNA/config/samples_metadata.txt",
                           header = TRUE)
# read.table(snakemake@input$meta, header = TRUE)
# Add library information
sample_table$libType <- "paired-end" 
# Simplify sample table
sample_table$rowname_IDs <- paste(sample_table$mutant, sample_table$run, "read_counts_sorted", sep = "/")
sample_table <- column_to_rownames(sample_table, var = "rowname_IDs")
sample_table_dex <- sample_table[,c("genotype", "libType")]
sample_table_dex <- dplyr::rename(sample_table_dex, condition = genotype)

# P635L
# 3 conditions
# Filter
P635L_count_files <- count_files[grep(pattern = "P635L", count_files)]
P635L_sample_table_dex <- sample_table_dex[grep(pattern = "P635L", row.names(sample_table_dex)),]
# Create DEXSeq object
P635L_dxd <- DEXSeqDataSetFromHTSeq(
  P635L_count_files,
  sampleData = P635L_sample_table_dex,
  design= ~ sample + exon + condition:exon,
  flattenedfile = flattened_file)
# the first 30 (we have 30 samples) corresponding to the number of reads mapping 
# to out exonic regions and the last seven correspond to the sum of the counts
# mapping to the rest of the exons from the same gene on each sample.
split(seq_len(ncol(P635L_dxd)), colData(P635L_dxd)$exon)
# rows are labelled with gene IDs , followed by a colon and the counting bin number
head(featureCounts(P635L_dxd), 5)
# details on the counting bins
head(rowRanges(P635L_dxd), 3)
# Design table
sampleAnnotation(P635L_dxd)

# Normalization
# DEXSeq uses the same method as DESeq2, which is provided in the function estimateSizeFactors().
P635L_dxd <- estimateSizeFactors(P635L_dxd)

# Dispersion estimation
# To test for differential exon usage, we need to estimate the variability of the data.
# The goal of this step is to distinguish technical and biological variation (noise) 
# from real effects on exon usage due to the different conditions
# TO NOTE: samples with the same experimental condition are considered as replicates 
P635L_dxd <- estimateDispersions(P635L_dxd)
plotDispEsts(P635L_dxd)

# Differential exon usage
# For each gene, DEXSeq fits a generalized linear model with the formula 
# ~sample + exon + condition:exon and compare it to the smaller model (the null model) ~ sample + exon
P635L_dxd <- testForDEU(P635L_dxd)

# From the coefficients of the fitted model, it is possible to distinguish overall 
# gene expression effects, that alter the counts from all the exons, 
# from exon usage effects, that are captured by the interaction term condition:exon 
# and that affect the each exon individually.
P635L_dxd <- estimateExonFoldChanges(P635L_dxd, fitExpToVar="condition")

# Results
P635L_dxd_res <- DEXSeqResults(P635L_dxd)
mcols(P635L_dxd_res)$description
# how many exonic regions are significant with a false discovery rate of 10%:
table(P635L_dxd_res$padj < 0.1)
# how many genes are affected
table(tapply(P635L_dxd_res$padj < 0.1, P635L_dxd_res$groupID, any))
# MA plot
# logarithm of fold change versus average normalized count per exon and marks 
# by red colour the exons which are considered significant
# --> how the power to detect differential exon usage depends on the number of reads that map to an exon
plotMA(P635L_dxd_res, cex=0.8)

# Visualization
P635L_dxd_res[which(P635L_dxd_res$padj < 0.1),]
plotDEXSeq(P635L_dxd_res, "ENSMUSG00000020359.14", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2)
# Ensembl info:
# http://www.ensembl.org/Mus_musculus/Transcript/Summary?g=ENSMUSG00000020359;r=11:51490094-51491894;t=ENSMUST00000169995

# hom vs WT
# Filter
P635L_sample_table_dex_wt_hom <- dplyr::filter(P635L_sample_table_dex, condition != "het")
idx_P635L_count_files_wt_hom <- P635L_count_files %in% paste0("./results/STAR/", rownames(P635L_sample_table_dex_wt_hom), ".txt")
P635L_count_files_wt_hom <- P635L_count_files[idx_P635L_count_files_wt_hom]
# Create DEXSeq object
P635L_dxd_wt_hom <- DEXSeqDataSetFromHTSeq(
  P635L_count_files_wt_hom,
  sampleData = P635L_sample_table_dex_wt_hom,
  design= ~ sample + exon + condition:exon,
  flattenedfile = flattened_file)
split(seq_len(ncol(P635L_dxd_wt_hom)), colData(P635L_dxd_wt_hom)$exon)

# Normalization
P635L_dxd_wt_hom <- estimateSizeFactors(P635L_dxd_wt_hom)

# Dispersion estimation
P635L_dxd_wt_hom <- estimateDispersions(P635L_dxd_wt_hom)
plotDispEsts(P635L_dxd_wt_hom)

# Differential exon usage
P635L_dxd_wt_hom <- testForDEU(P635L_dxd_wt_hom)

# Estimate exon fold change
P635L_dxd_wt_hom <- estimateExonFoldChanges(P635L_dxd_wt_hom, fitExpToVar="condition")

# Results
P635L_dxd_res_wt_hom <- DEXSeqResults(P635L_dxd_wt_hom)
mcols(P635L_dxd_res_wt_hom)$description
# how many exonic regions are significant with a false discovery rate of 10%:
table(P635L_dxd_res_wt_hom$padj < 0.1)
# how many genes are affected
table(tapply(P635L_dxd_res_wt_hom$padj < 0.1, P635L_dxd_res_wt_hom$groupID, any))
# MA plot
plotMA(P635L_dxd_res_wt_hom, cex=0.8)


# R636Q
# 3 conditions
# Filter
R636Q_count_files <- count_files[grep(pattern = "R636Q", count_files)]
R636Q_sample_table_dex <- sample_table_dex[grep(pattern = "R636Q", row.names(sample_table_dex)),]
# Create DEXSeq object
R636Q_dxd <- DEXSeqDataSetFromHTSeq(
  R636Q_count_files,
  sampleData = R636Q_sample_table_dex,
  design= ~ sample + exon + condition:exon,
  flattenedfile = flattened_file)
split(seq_len(ncol(R636Q_dxd)), colData(R636Q_dxd)$exon)

# Normalization
R636Q_dxd <- estimateSizeFactors(R636Q_dxd)

# Dispersion estimation
R636Q_dxd <- estimateDispersions(R636Q_dxd)
plotDispEsts(R636Q_dxd)

# Differential exon usage
R636Q_dxd <- testForDEU(R636Q_dxd)

# Exon Fold Changes
R636Q_dxd <- estimateExonFoldChanges(R636Q_dxd, fitExpToVar="condition")

# Results
R636Q_dxd_res <- DEXSeqResults(R636Q_dxd)
mcols(R636Q_dxd_res)$description
# how many exonic regions are significant with a false discovery rate of 10%:
table(R636Q_dxd_res$padj < 0.1)
# how many genes are affected
table(tapply(R636Q_dxd_res$padj < 0.1, R636Q_dxd_res$groupID, any))
# MA plot
plotMA(R636Q_dxd_res, cex=0.8)

# Visualization
R636Q_dxd_res[which(R636Q_dxd_res$padj < 0.1),]
plotDEXSeq(R636Q_dxd_res, "ENSMUSG00000020359.14", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2)
# Ensembl info:
# http://www.ensembl.org/Mus_musculus/Transcript/Summary?g=ENSMUSG00000020359;r=11:51490094-51491894;t=ENSMUST00000169995

# hom vs WT
# Filter
R636Q_sample_table_dex_wt_hom <- dplyr::filter(R636Q_sample_table_dex, condition != "het")
idx_R636Q_count_files_wt_hom <- R636Q_count_files %in% paste0("./results/STAR/", rownames(R636Q_sample_table_dex_wt_hom), ".txt")
R636Q_count_files_wt_hom <- R636Q_count_files[idx_R636Q_count_files_wt_hom]
# Create DEXSeq object
R636Q_dxd_wt_hom <- DEXSeqDataSetFromHTSeq(
  R636Q_count_files_wt_hom,
  sampleData = R636Q_sample_table_dex_wt_hom,
  design= ~ sample + exon + condition:exon,
  flattenedfile = flattened_file)
split(seq_len(ncol(R636Q_dxd_wt_hom)), colData(R636Q_dxd_wt_hom)$exon)

# Normalization
R636Q_dxd_wt_hom <- estimateSizeFactors(R636Q_dxd_wt_hom)

# Dispersion estimation
R636Q_dxd_wt_hom <- estimateDispersions(R636Q_dxd_wt_hom)
plotDispEsts(R636Q_dxd_wt_hom)

# Differential exon usage
R636Q_dxd_wt_hom <- testForDEU(R636Q_dxd_wt_hom)

# Estimate exon fold change
R636Q_dxd_wt_hom <- estimateExonFoldChanges(R636Q_dxd_wt_hom, fitExpToVar="condition")

# Results
R636Q_dxd_res_wt_hom <- DEXSeqResults(R636Q_dxd_wt_hom)
mcols(R636Q_dxd_res_wt_hom)$description
# how many exonic regions are significant with a false discovery rate of 10%:
table(R636Q_dxd_res_wt_hom$padj < 0.1)
# how many genes are affected
table(tapply(R636Q_dxd_res_wt_hom$padj < 0.1, R636Q_dxd_res_wt_hom$groupID, any))
# MA plot
plotMA(R636Q_dxd_res_wt_hom, cex=0.8)


# Read gtf file
# gen_vM29_annot <- read.delim("resources/gencode.vM29.annotation.gtf", header=F, comment.char="#")
gen_vM29_annot <- import("resources/gencode.vM29.annotation.gtf")
# Extract info
gen_vM29_annot_gene <- gen_vM29_annot[,c("gene_id", "gene_name")]
gen_vM29_annot_gene_uniq <- unique(as.data.frame(gen_vM29_annot_gene)[,c("gene_id", "gene_name")])
colnames(gen_vM29_annot_gene_uniq)[1] <- c("groupID")
# Overlap results
# P635L WT vs HOM
# P635L_dxd_res_wt_hom
# P635L_dxd_res_wt_hom$groupID
df_P635L_dxd_res_wt_hom_gene_name <- left_join(as.data.frame(P635L_dxd_res_wt_hom), 
                                               gen_vM29_annot_gene_uniq)
# Plot
# Ttn
P635L_idx_Ttn <- which(df_P635L_dxd_res_wt_hom_gene_name$gene_name == "Ttn")
df_P635L_dxd_res_wt_hom_gene_name[P635L_idx_Ttn,]
P635L_Ttn_ensembl <- unique(P635L_dxd_res_wt_hom[P635L_idx_Ttn,"groupID"])
plotDEXSeq(P635L_dxd_res_wt_hom, 
           P635L_Ttn_ensembl, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2)
P635L_dxd_res_wt_hom[P635L_idx_Ttn,]
as.data.frame(P635L_dxd_res_wt_hom[P635L_idx_Ttn,c("featureID", "pvalue", "padj")])
# Ryr2
P635L_idx_Ryr2 <- which(df_P635L_dxd_res_wt_hom_gene_name$gene_name == "Ryr2")
df_P635L_dxd_res_wt_hom_gene_name[P635L_idx_Ryr2,]
P635L_Ryr2_ensembl <- unique(P635L_dxd_res_wt_hom[P635L_idx_Ryr2,"groupID"])
plotDEXSeq(P635L_dxd_res_wt_hom, 
           P635L_Ryr2_ensembl, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2)
P635L_dxd_res_wt_hom[P635L_idx_Ryr2,]
as.data.frame(P635L_dxd_res_wt_hom[P635L_idx_Ryr2,c("featureID", "pvalue", "padj")])
# Ldb3
P635L_idx_Ldb3 <- which(df_P635L_dxd_res_wt_hom_gene_name$gene_name == "Ldb3")
df_P635L_dxd_res_wt_hom_gene_name[P635L_idx_Ldb3,]
P635L_Ldb3_ensembl <- unique(P635L_dxd_res_wt_hom[P635L_idx_Ldb3,"groupID"])
plotDEXSeq(P635L_dxd_res_wt_hom, 
           P635L_Ldb3_ensembl, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2)
P635L_dxd_res_wt_hom[P635L_idx_Ldb3,]
as.data.frame(P635L_dxd_res_wt_hom[P635L_idx_Ldb3,c("featureID", "pvalue", "padj")])
# Camk2d
P635L_idx_Camk2d <- which(df_P635L_dxd_res_wt_hom_gene_name$gene_name == "Camk2d")
df_P635L_dxd_res_wt_hom_gene_name[P635L_idx_Camk2d,]
P635L_Camk2d_ensembl <- unique(P635L_dxd_res_wt_hom[P635L_idx_Camk2d,"groupID"])
plotDEXSeq(P635L_dxd_res_wt_hom, 
           P635L_Camk2d_ensembl, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2)
P635L_dxd_res_wt_hom[P635L_idx_Camk2d,]
as.data.frame(P635L_dxd_res_wt_hom[P635L_idx_Camk2d,c("featureID", "pvalue", "padj")])


# R636Q WT vs HOM
# R636Q_dxd_res_wt_hom
# R636Q_dxd_res_wt_hom$groupID
df_R636Q_dxd_res_wt_hom_gene_name <- left_join(as.data.frame(R636Q_dxd_res_wt_hom), 
                                               gen_vM29_annot_gene_uniq)
# Plot
# Ttn
R636Q_idx_Ttn <- which(df_R636Q_dxd_res_wt_hom_gene_name$gene_name == "Ttn")
df_R636Q_dxd_res_wt_hom_gene_name[R636Q_idx_Ttn,]
R636Q_Ttn_ensembl <- unique(R636Q_dxd_res_wt_hom[R636Q_idx_Ttn,"groupID"])
plotDEXSeq(R636Q_dxd_res_wt_hom, 
           R636Q_Ttn_ensembl, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2, norCounts=TRUE, splicing=TRUE)
R636Q_dxd_res_wt_hom[R636Q_idx_Ttn,]
as.data.frame(R636Q_dxd_res_wt_hom[R636Q_idx_Ttn,c("featureID", "pvalue", "padj")])
# Ryr2
R636Q_idx_Ryr2 <- which(df_R636Q_dxd_res_wt_hom_gene_name$gene_name == "Ryr2")
df_R636Q_dxd_res_wt_hom_gene_name[R636Q_idx_Ryr2,]
R636Q_Ryr2_ensembl <- unique(R636Q_dxd_res_wt_hom[R636Q_idx_Ryr2,"groupID"])
plotDEXSeq(R636Q_dxd_res_wt_hom, 
           R636Q_Ryr2_ensembl, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2)
R636Q_dxd_res_wt_hom[R636Q_idx_Ryr2,]
as.data.frame(R636Q_dxd_res_wt_hom[R636Q_idx_Ryr2,c("featureID", "pvalue", "padj")])
# Ldb3
R636Q_idx_Ldb3 <- which(df_R636Q_dxd_res_wt_hom_gene_name$gene_name == "Ldb3")
df_R636Q_dxd_res_wt_hom_gene_name[R636Q_idx_Ldb3,]
R636Q_Ldb3_ensembl <- unique(R636Q_dxd_res_wt_hom[R636Q_idx_Ldb3,"groupID"])
plotDEXSeq(R636Q_dxd_res_wt_hom, 
           R636Q_Ldb3_ensembl, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2,norCounts=TRUE, splicing=TRUE)
R636Q_dxd_res_wt_hom[R636Q_idx_Ldb3,]
as.data.frame(R636Q_dxd_res_wt_hom[R636Q_idx_Ldb3,c("featureID", "pvalue", "padj")])
# Camk2d
R636Q_idx_Camk2d <- which(df_R636Q_dxd_res_wt_hom_gene_name$gene_name == "Camk2d")
df_R636Q_dxd_res_wt_hom_gene_name[R636Q_idx_Camk2d,]
R636Q_Camk2d_ensembl <- unique(R636Q_dxd_res_wt_hom[R636Q_idx_Camk2d,"groupID"])
plotDEXSeq(R636Q_dxd_res_wt_hom, 
           R636Q_Camk2d_ensembl, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2)
plotDEXSeq(R636Q_dxd_res_wt_hom, 
           R636Q_Camk2d_ensembl, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2, norCounts=TRUE)
plotDEXSeq(R636Q_dxd_res_wt_hom, 
           R636Q_Camk2d_ensembl, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2,
           norCounts=TRUE, splicing=TRUE)

R636Q_dxd_res_wt_hom[R636Q_idx_Camk2d,]
as.data.frame(R636Q_dxd_res_wt_hom[R636Q_idx_Camk2d,c("featureID", "pvalue", "padj")])


# NAs
# GLM fit that calculates the effect sizes in the estimateExonFoldChanges fails to converge
# --> did not see an error when running the code
# insufficient library size / coverage
hist(R636Q_dxd_res_wt_hom$exonBaseMean, breaks=100)
hist( log2(1 + R636Q_dxd_res_wt_hom$exonBaseMean), breaks=100)
plotMA(R636Q_dxd_res_wt_hom) 
# scatter plot of log2 fold changes (on the y-axis) versus the mean of normalized counts (on the x-axis).

hist(R636Q_dxd_res_wt_hom$pvalue)
table(R636Q_dxd_res_wt_hom$pvalue < 0.05)

# Plot the lowest p-values
min(R636Q_dxd_res_wt_hom$pvalue, na.rm = TRUE)
R636Q_dxd_res_wt_hom[order(R636Q_dxd_res_wt_hom$pvalue)[1:3],]
R636Q_dxd_res_wt_hom[order(R636Q_dxd_res_wt_hom$pvalue)[1:3],"groupID"]
idx_R636Q_wt_hom_smallest <- order(R636Q_dxd_res_wt_hom$pvalue)[1:3]
df_R636Q_dxd_res_wt_hom_gene_name[idx_R636Q_wt_hom_smallest, "gene_name"]
plotDEXSeq(R636Q_dxd_res_wt_hom, 
           R636Q_dxd_res_wt_hom[idx_R636Q_wt_hom_smallest[1], "groupID"], 
           legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2, norCounts=TRUE, splicing=TRUE)
#the script might come across overlapping genes. If two genes on the same strand are 
#found with an exon of the first gene overlapping with an exon of the second gene, 
#the script’s default behavior is to combine the genes into a single “aggregate gene” 
#which is subsequently referred to with the IDs of the individual genes, joined by a plus (‘+’) sign.
plotDEXSeq(R636Q_dxd_res_wt_hom, 
           R636Q_dxd_res_wt_hom[idx_R636Q_wt_hom_smallest[2], "groupID"], 
           legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2, norCounts=TRUE, splicing=TRUE)
plotDEXSeq(R636Q_dxd_res_wt_hom, 
           R636Q_dxd_res_wt_hom[idx_R636Q_wt_hom_smallest[3], "groupID"], 
           legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2, norCounts=TRUE, splicing=TRUE)
# Filter no exon?
count_files[16]
test <- read.table(count_files[16])
head(test)
unique(test$V2)
table(test$V2)
test[1:9,]
R636Q_dxd_res_wt_hom[1:9,]

head( featureCounts(R636Q_dxd), 50)
# Use the parameter 'minCount' in estimateDispersions for this. The
# default, 10, is probably reasonable, but maybe a bit on the low end.
# Every exon with less than 10 counts (summed across all sample) is
# omitted from the tests.

# check number of reads
R636Q_Camk2d_ensembl
grep(pattern = paste0('"',R636Q_Camk2d_ensembl, ':*"'), test$V1)

# Remove zero counts genes
assay(R636Q_dxd_wt_hom)
colData(R636Q_dxd_wt_hom)
idx_R636Q_dxd_wt_hom_diffzero <- which(rowSums(assay(R636Q_dxd_wt_hom)) != 0)
R636Q_dxd_wt_hom_diffzero <- R636Q_dxd_wt_hom[idx_R636Q_dxd_wt_hom_diffzero,]
# Normalization
R636Q_dxd_wt_hom_diffzero <- estimateSizeFactors(R636Q_dxd_wt_hom_diffzero)
# Dispersion estimation
R636Q_dxd_wt_hom_diffzero <- estimateDispersions(R636Q_dxd_wt_hom_diffzero)
plotDispEsts(R636Q_dxd_wt_hom_diffzero)
# Differential exon usage
R636Q_dxd_wt_hom_diffzero <- testForDEU(R636Q_dxd_wt_hom_diffzero)
# Estimate exon fold change
R636Q_dxd_wt_hom_diffzero <- estimateExonFoldChanges(R636Q_dxd_wt_hom_diffzero, 
                                                     fitExpToVar="condition")
# Results
R636Q_dxd_res_wt_hom_diffzero <- DEXSeqResults(R636Q_dxd_wt_hom_diffzero)
mcols(R636Q_dxd_res_wt_hom_diffzero)$description
# how many exonic regions are significant with a false discovery rate of 10%:
table(R636Q_dxd_res_wt_hom_diffzero$padj < 0.1)
# how many genes are affected
table(tapply(R636Q_dxd_res_wt_hom_diffzero$padj < 0.1, R636Q_dxd_res_wt_hom_diffzero$groupID, any))
# MA plot
plotMA(R636Q_dxd_res_wt_hom_diffzero, cex=0.8)

# Filter and apply pval correction only for samples with more than 20 reads
R636Q_dxd_res_wt_hom$pvalue[rowSums( R636Q_dxd_res_wt_hom$countData ) < 20] <- NA
R636Q_dxd_res_wt_hom$padj <- p.adjust(R636Q_dxd_res_wt_hom$pvalue, method="BH")

# Number of reads
colSums(counts(R636Q_dxd_res_wt_hom))
counts(R636Q_dxd_res_wt_hom)
# R636Q/22s000641/reads_counts_clean 9407
# /g/steinmetz/project/dcm_lgreads/mouse_bulkRNA/results/STAR/R636Q/22s000641
# awk -F'\t' '{sum+=$2;} END{print sum;}' reads_counts.txt
# 336001

# https://www.biostars.org/p/405934/

# DESeq number of reads
library(DESeq2)
# ?DESeqDataSetFromHTSeq
library(Rsubread)
# featureCounts() from the Rsubread package
featurecounts_22s000641 <- featureCounts(file = "./results/STAR/R636Q/22s000641/Aligned.sortedByCoord.out.bam",
              annot.ext = "resources/gencode.vM29.annotation.gtf",
              isGTFAnnotationFile = TRUE,
              strandSpecific = 1,
              isPairedEnd = TRUE)
colSums(featurecounts_22s000641$counts)
head(DEXSeq::featureCounts(R636Q_dxd_wt_hom), 5)
sum(DEXSeq::featureCounts(R636Q_dxd_wt_hom)[,"R636Q/22s000641/reads_counts_clean"])
# sort bam file

################################################################################
# Thread to import from Kallisto to DEXseq:
## https://support.bioconductor.org/p/113693/
# Same code as for the DESeq analysis #
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
                    tx2gene = unique(tx_g_ids), 
                    txOut=TRUE) # to get transcripts level counts matrix

# Counts data extraction
counts_all_mice <- t(all_txi$counts)
counts_all_mice_wnames <- rownames_to_column(as.data.frame(counts_all_mice), 
                                             var = "run")
counts_all_mice_wmeta <- left_join(counts_all_mice_wnames, 
                                   samples_metadata)

# P635L
# Filter 
# Counts
idx_P635L <- which(counts_all_mice_wmeta$mutant == "P635L")
counts_P635L_mice <- counts_all_mice[idx_P635L,]
t_counts_P635L_mice <- t(counts_P635L_mice)
# Match transcript IDs to gene IDs
# dim(t_counts_P635L_mice)
# row_names_ensmus <- sapply(row.names(t_counts_P635L_mice), FUN = function(x){
#   unlist(strsplit(x, "[.]"))[1]
#   })
# length(unique(row_names_ensmus))
group_id <- left_join(data.frame("TXNAME" = row.names(t_counts_P635L_mice)),
          tx_g_ids)
# Filter transcript without gene names
idx_notna <- which(group_id$GENEID != "NA")
# Metadata
samples_metadata_P635L <- dplyr::filter(samples_metadata,
                                        mutant == "P635L")
# DEXSeq
dexseq_P635L <- DEXSeqDataSet(countData = round(t_counts_P635L_mice)[idx_notna,],
                              sampleData = samples_metadata_P635L[,c("run", "genotype")],
                              design = ~ run + genotype + exon:genotype, 
                              featureID = row.names(t_counts_P635L_mice[idx_notna,]),
                              groupID = group_id$GENEID[idx_notna])
colData(dexseq_P635L)
head(counts(dexseq_P635L), 5)
head(featureCounts(dexseq_P635L), 5)
head(rowRanges(dexseq_P635L), 3)
sampleAnnotation(dexseq_P635L)
# gff_file <- read.delim("results/Mus_musculus.GRCm38.96.DEXSeq.chr.gff", 
#                        header=F, 
#                        comment.char="#") 
# Sample normalization
dexseq_P635L <- estimateSizeFactors(dexseq_P635L)
# Dispersion estimation
dexseq_P635L <- estimateDispersions(dexseq_P635L)
plotDispEsts(dexseq_P635L)
# Testing for differential exon usage
dexseq_P635L <- testForDEU(dexseq_P635L)
# estimate relative exon usage fold changes
dexseq_P635L <- estimateExonFoldChanges(dexseq_P635L, fitExpToVar="genotype")
# summarize the results
dexseq_P635L_res <- DEXSeqResults(dexseq_P635L)
mcols(dexseq_P635L_res)$description
# how many exonic regions are significant with a false discovery rate of 10%:
table(dexseq_P635L_res$padj < 0.1)
# how many genes are affected
table(tapply(dexseq_P635L_res$padj < 0.1, dexseq_P635L_res$groupID, any))
#  detect differential exon usage depends on the number of reads that map to an exon, 
# a so-called MA plot is useful, which plots the logarithm of fold change versus 
# average normalized count per exon and marks by red colour the exons which are considered 
# significant
plotMA(dexseq_P635L_res, cex=0.8)



