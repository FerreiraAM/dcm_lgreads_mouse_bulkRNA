# DEXSeq perform DEU test

# Load libraries
library(DEXSeq)

# Mutation
mut <- snakemake@wildcards$mutation

# Load data
load(file = paste0("results/DEXSeq/robject_", {mut}, "_dexseq_estdisp_WTvsHOM.rda"))
# Differential exon usage
dxd_wt_hom <- testForDEU(dxd_wt_hom)

# Estimate exon fold change
dxd_wt_hom <- estimateExonFoldChanges(dxd_wt_hom, fitExpToVar="condition")

# Results
dxd_res_wt_hom <- DEXSeqResults(dxd_wt_hom)
# Save object
save(dxd_res_wt_hom, file = paste0("results/DEXSeq/robject_", {mut}, "_dexseq_DEUres_WTvsHOM.rda"))