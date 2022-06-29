# DEXSeq analysis
# Thread to import from Kallisto to DEXseq:
## https://support.bioconductor.org/p/113693/


# Load libraries
library(DEXSeq)
library(tximport)

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
# Metadata
samples_metadata_P635L <- dplyr::filter(samples_metadata,
                                        mutant == "P635L")
# DEXSeq
DEXSeqDataSet(countData = round(t_counts_P635L_mice),
              sampleData = samples_metadata_P635L[,c("run", "genotype")],
              design = ~ run + genotype + exon:genotype, 
              featureID = row.names(t_counts_P635L_mice),
              groupID = row.names(t_counts_P635L_mice))
