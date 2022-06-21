# DEBUG SCRIPT
# save.image("deseq_analysis.rda")
# stop()

# Load libraries:
library(tidyverse) 
library(DESeq2)
library(tximport)
library(rtracklayer)
library(ggplot2)
library(ggfortify)

# Data analysis 

# Gene ID data frame
tx_g_ids <- read.table("/g/steinmetz/project/dcm_lgreads/mouse_bulkRNA/resources/mus_musculus/transcripts_to_genes.txt")
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
          file = "/g/steinmetz/project/dcm_lgreads/mouse_bulkRNA/results/total_nb_reads.csv", 
          col.names = TRUE, row.names = FALSE)

# PCA
# Prepare data for PCA
counts_all_mice <- t(all_txi$counts)
counts_all_mice_wnames <- rownames_to_column(as.data.frame(counts_all_mice), var = "run")
counts_all_mice_wmeta <- left_join(counts_all_mice_wnames, samples_metadata)
# Perform PCA
pca_all_mice <- prcomp(counts_all_mice)
plot_pca <- autoplot(pca_all_mice, data = counts_all_mice_wmeta, colour = "mutant") +
  theme_bw() +
  ggtitle("PCA with all mice and all genes")
ggsave("/g/steinmetz/project/dcm_lgreads/mouse_bulkRNA/results/plot_pca_all_mice_genes.pdf", 
       plot = plot_pca)

# DESeq2
