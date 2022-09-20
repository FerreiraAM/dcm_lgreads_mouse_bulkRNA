# To visualize DEXSeq DEU results

# Load libraries
library(DEXSeq)
library(rtracklayer)
library(tidyverse)

# Mutation
# mut <- snakemake@wildcards$mutation

# Number of reads
sample_table <- read.table("/g/steinmetz/project/dcm_lgreads/mouse_bulkRNA/config/samples_metadata.txt",
                           header = TRUE)
# read.table(snakemake@input$meta, header = TRUE)
sample_table$samples <- paste(sample_table$mutant, sample_table$run, "read_counts_sorted_reverse_clean", sep = "/")
# Filter het
sample_table <- dplyr::filter(sample_table, genotype != "het")
# P635L DEXSeq object
load(file = paste0("results/DEXSeq/robject_P635L_dexseq_estdisp_WTvsHOM.rda"))
# P635L Read counts
P635L_counts <- colSums(DEXSeq::featureCounts(dxd_wt_hom))
P635L_total_readcounts_DEXSeq <- data.frame(samples=names(P635L_counts), 
                                            read_counts=P635L_counts, 
                                            row.names=NULL)
rm(dxd_wt_hom)
# R636Q DEXSeq object
load(file = paste0("results/DEXSeq/robject_R636Q_dexseq_estdisp_WTvsHOM.rda"))
# R636Q Read counts
R636Q_counts <- colSums(DEXSeq::featureCounts(dxd_wt_hom))
R636Q_total_readcounts_DEXSeq <- data.frame(samples=names(R636Q_counts), 
                                            read_counts=R636Q_counts, 
                                            row.names=NULL)
# Merge
P635L_R636Q_total_readcounts_DEXSeq <- rbind(P635L_total_readcounts_DEXSeq, 
                                       R636Q_total_readcounts_DEXSeq)
P635L_R636Q_total_readcounts_DEXSeq <- left_join(sample_table, 
                                                 P635L_R636Q_total_readcounts_DEXSeq)
write.table(P635L_R636Q_total_readcounts_DEXSeq[,c("mutant", "run", "genotype", "read_counts")], 
            file = paste0("results/DEXSeq/P635L_R636Q_total_readcounts_DEXSeq.txt"),
            col.names = TRUE,
            row.names = FALSE)

# Read gtf file
gen_vM29_annot <- import("resources/gencode.vM29.annotation.gtf")
# Extract info
gen_vM29_annot_gene <- gen_vM29_annot[,c("gene_id", "gene_name")]
gen_vM29_annot_gene_uniq <- unique(as.data.frame(gen_vM29_annot_gene)[,c("gene_id", "gene_name")])
colnames(gen_vM29_annot_gene_uniq)[1] <- c("groupID")

# Load DEU results
# P635L
name_P635L_dxd_res_wt_hom <- load(file = paste0("results/DEXSeq/robject_P635L_dexseq_DEUres_WTvsHOM.rda"))
# Get the object by its name
P635L_dxd_res_wt_hom <- get(name_P635L_dxd_res_wt_hom)
# Remove the old object
rm(dxd_res_wt_hom)
# R636Q
name_R636Q_dxd_res_wt_hom <- load(file = paste0("results/DEXSeq/robject_R636Q_dexseq_DEUres_WTvsHOM.rda"))
# Get the object by its name
R636Q_dxd_res_wt_hom <- get(name_R636Q_dxd_res_wt_hom)
# Remove the old object
rm(dxd_res_wt_hom)

# Results data
# how many exonic regions are significant with a false discovery rate of 10%:
table(P635L_dxd_res_wt_hom$padj < 0.1)
table(R636Q_dxd_res_wt_hom$padj < 0.1)
# how many genes are affected
table(tapply(P635L_dxd_res_wt_hom$padj < 0.1, P635L_dxd_res_wt_hom$groupID, any))
v_P635L_groupID_sig <- dplyr::filter(as_tibble(P635L_dxd_res_wt_hom), padj < 0.1) %>% 
  pull(groupID) %>% 
  unique
table(tapply(R636Q_dxd_res_wt_hom$padj < 0.1, R636Q_dxd_res_wt_hom$groupID, any))
v_R636Q_groupID_sig <- dplyr::filter(as_tibble(R636Q_dxd_res_wt_hom), padj < 0.1) %>% 
  pull(groupID) %>% 
  unique
# MA plot
# logarithm of fold change versus average normalized count per exon and marks 
# by red colour the exons which are considered significant
# --> how the power to detect differential exon usage depends on the number of reads that map to an exon
plotMA(P635L_dxd_res_wt_hom, cex=0.8)
plotMA(R636Q_dxd_res_wt_hom, cex=0.8)

# Gene of interests
df_P635L_dxd_res_wt_hom_gene_name <- left_join(as.data.frame(P635L_dxd_res_wt_hom), 
                                               gen_vM29_annot_gene_uniq)
df_P635L_dxd_res_wt_hom_gene_name_sig <- dplyr::filter(df_P635L_dxd_res_wt_hom_gene_name, 
                                                       padj < 0.1) %>% 
  dplyr::select(groupID, featureID, pvalue, padj, log2fold_wt_hom, log2fold_wt_hom, gene_name)
write.table(df_P635L_dxd_res_wt_hom_gene_name_sig, 
            file = "results/DEXSeq/table_P635L_dxd_res_wt_hom_gene_name_sig.txt",
            col.names = TRUE,
            row.names = FALSE)
df_R636Q_dxd_res_wt_hom_gene_name <- left_join(as.data.frame(R636Q_dxd_res_wt_hom), 
                                               gen_vM29_annot_gene_uniq)
# Ttn
P635L_idx_Ttn <- which(df_P635L_dxd_res_wt_hom_gene_name$gene_name == "Ttn")
R636Q_idx_Ttn <- which(df_R636Q_dxd_res_wt_hom_gene_name$gene_name == "Ttn")
Ttn_ensembl <- unique(R636Q_dxd_res_wt_hom[R636Q_idx_Ttn,"groupID"])
# P635L plot
plotDEXSeq(P635L_dxd_res_wt_hom, 
           Ttn_ensembl, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2, 
           norCounts=TRUE, splicing=TRUE)
P635L_dxd_res_wt_hom[P635L_idx_Ttn,]
as.data.frame(P635L_dxd_res_wt_hom[P635L_idx_Ttn,c("featureID", "pvalue", "padj")])
# R636Q plot
plotDEXSeq(R636Q_dxd_res_wt_hom, 
           Ttn_ensembl, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2, 
           norCounts=TRUE, splicing=TRUE)
R636Q_dxd_res_wt_hom[R636Q_idx_Ttn,]
as.data.frame(R636Q_dxd_res_wt_hom[R636Q_idx_Ttn,c("featureID", "pvalue", "padj")])
# Ryr2
P635L_idx_Ryr2 <- which(df_P635L_dxd_res_wt_hom_gene_name$gene_name == "Ryr2")
R636Q_idx_Ryr2 <- which(df_R636Q_dxd_res_wt_hom_gene_name$gene_name == "Ryr2")
Ryr2_ensembl <- unique(R636Q_dxd_res_wt_hom[R636Q_idx_Ryr2,"groupID"])
# P635L plot
plotDEXSeq(P635L_dxd_res_wt_hom,
           Ryr2_ensembl, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2,
           norCounts=TRUE, splicing=TRUE)
# glm fit fails and the coefficients are not estimated
P635L_dxd_res_wt_hom[P635L_idx_Ryr2,]
as.data.frame(P635L_dxd_res_wt_hom[P635L_idx_Ryr2,c("featureID", "pvalue", "padj")])
# R636Q plot
plotDEXSeq(R636Q_dxd_res_wt_hom, 
           Ryr2_ensembl, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2, 
           norCounts=TRUE, splicing=TRUE)
R636Q_dxd_res_wt_hom[R636Q_idx_Ryr2,]
as.data.frame(R636Q_dxd_res_wt_hom[R636Q_idx_Ryr2,c("featureID", "pvalue", "padj")])
# Ldb3
P635L_idx_Ldb3 <- which(df_P635L_dxd_res_wt_hom_gene_name$gene_name == "Ldb3")
R636Q_idx_Ldb3 <- which(df_R636Q_dxd_res_wt_hom_gene_name$gene_name == "Ldb3")
Ldb3_ensembl <- unique(R636Q_dxd_res_wt_hom[R636Q_idx_Ldb3,"groupID"])
# P635L plot
plotDEXSeq(P635L_dxd_res_wt_hom, 
           Ldb3_ensembl, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2, 
           norCounts=TRUE, splicing=TRUE)
P635L_dxd_res_wt_hom[P635L_idx_Ldb3,]
as.data.frame(P635L_dxd_res_wt_hom[P635L_idx_Ldb3,c("featureID", "pvalue", "padj")])
# R636Q plot
plotDEXSeq(R636Q_dxd_res_wt_hom, 
           Ldb3_ensembl, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2,
           norCounts=TRUE, splicing=TRUE)
R636Q_dxd_res_wt_hom[R636Q_idx_Ldb3,]
as.data.frame(R636Q_dxd_res_wt_hom[R636Q_idx_Ldb3,c("featureID", "pvalue", "padj")])
# Camk2d
P635L_idx_Camk2d <- which(df_P635L_dxd_res_wt_hom_gene_name$gene_name == "Camk2d")
R636Q_idx_Camk2d <- which(df_R636Q_dxd_res_wt_hom_gene_name$gene_name == "Camk2d")
Camk2d_ensembl <- unique(R636Q_dxd_res_wt_hom[R636Q_idx_Camk2d,"groupID"])
# P635L plot
plotDEXSeq(P635L_dxd_res_wt_hom, 
           Camk2d_ensembl, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2, 
           norCounts=TRUE, splicing=TRUE)
P635L_dxd_res_wt_hom[P635L_idx_Camk2d,]
as.data.frame(P635L_dxd_res_wt_hom[P635L_idx_Camk2d,c("featureID", "pvalue", "padj")])
# R636Q plot
plotDEXSeq(R636Q_dxd_res_wt_hom, 
           Camk2d_ensembl, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2,
           norCounts=TRUE, splicing=TRUE)
R636Q_dxd_res_wt_hom[R636Q_idx_Camk2d,]
as.data.frame(R636Q_dxd_res_wt_hom[R636Q_idx_Camk2d,c("featureID", "pvalue", "padj")])


dplyr::filter(df_R636Q_dxd_res_wt_hom_gene_name, padj < 0.1)

############## DEBUG R636Q ##############
# Some sample labels might have been switched
library(ggrepel)
# Ldb3
check_Ldb3_values <- dplyr::filter(df_R636Q_dxd_res_wt_hom_gene_name, 
                                   gene_name == "Ldb3")
dim(dplyr::select(check_Ldb3_values, contains("countData")))

plotDEXSeq(R636Q_dxd_res_wt_hom, 
           Ldb3_ensembl, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2,
           norCounts=TRUE, expression = FALSE)
# Long format
check_Ldb3_counts_exons_samp <- dplyr::select(check_Ldb3_values, 
                                              featureID, contains("countData"))
check_Ldb3_counts_exons_samp_lg <- pivot_longer(check_Ldb3_counts_exons_samp, 
             cols = contains("countData")) %>% 
  mutate(run = sub(".+R636Q\\.(.+)\\.read_counts.+", "\\1", name))
check_Ldb3_counts_exons_samp_lg <- check_Ldb3_counts_exons_samp_lg %>% 
  left_join(sample_table) %>% 
  unite(sample_genotype, run, genotype, sep = "_", remove = FALSE)
# Plot
ggplot(check_Ldb3_counts_exons_samp_lg, 
       aes(x = featureID, 
           y = value, 
           color = as.factor(genotype), 
           shape = as.factor(sample_genotype))) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_shape_manual(values=seq(0,15)) +
  ggtitle("R636Q - Ldb3 gene")
# Select exon E16 to E20
dplyr::filter(check_Ldb3_counts_exons_samp_lg, featureID %in% paste0("E0", 16:20)) %>% 
  ggplot(., aes(x = featureID, 
                y = value, 
                color = as.factor(genotype), 
                shape = as.factor(sample_genotype),
                label = sample_genotype)) +
  geom_point(cex = 3) +
  geom_text(hjust=0, vjust=0.2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_shape_manual(values=seq(0,15)) +
  ggtitle("R636Q - Ldb3 gene")
# Ttn
check_Ttn_values <- dplyr::filter(df_R636Q_dxd_res_wt_hom_gene_name, 
                                   gene_name == "Ttn")
dim(dplyr::select(check_Ttn_values, contains("countData")))
plotDEXSeq(R636Q_dxd_res_wt_hom, 
           Ttn_ensembl, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2,
           norCounts=TRUE, expression = FALSE)
check_Ttn_counts_exons_samp <- dplyr::select(check_Ttn_values, 
                                              featureID, contains("countData"))
check_Ttn_counts_exons_samp_lg <- pivot_longer(check_Ttn_counts_exons_samp, 
                                                cols = contains("countData")) %>% 
  mutate(run = sub(".+R636Q\\.(.+)\\.read_counts.+", "\\1", name))
check_Ttn_counts_exons_samp_lg <- check_Ttn_counts_exons_samp_lg %>% 
  left_join(sample_table) %>% 
  unite(sample_genotype, run, genotype, sep = "_", remove = FALSE)
# Plot
ggplot(check_Ttn_counts_exons_samp_lg, 
       aes(x = featureID, 
           y = value, 
           color = as.factor(genotype), 
           shape = as.factor(sample_genotype))) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_shape_manual(values=seq(0,15)) +
  ggtitle("R636Q - Ttn gene")
# Select exon E16 to E20
dplyr::filter(check_Ttn_counts_exons_samp_lg, featureID %in% paste0("E", 290:330)) %>% 
  ggplot(., aes(x = featureID, 
                y = value, 
                color = as.factor(genotype), 
                shape = as.factor(sample_genotype))) +
  geom_point(cex = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_shape_manual(values=seq(0,15)) +
  ggtitle("R636Q - Ttn gene")
dplyr::filter(check_Ttn_counts_exons_samp_lg, featureID %in% paste0("E", 300:310)) %>% 
  ggplot(., aes(x = featureID, 
                y = value, 
                color = as.factor(genotype), 
                shape = as.factor(sample_genotype),
                label = sample_genotype)) +
  geom_point(cex = 3) +
  geom_text(hjust=0, vjust=0.2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_shape_manual(values=seq(0,15)) +
  ggtitle("R636Q - Ttn gene")
