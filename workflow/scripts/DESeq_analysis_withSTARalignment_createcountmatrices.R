# Create count matrices from STAR alignment for DESeq2 analysis 
# save.image("alignment.rda")
# stop()

# Load libraries
library(Rsubread)

# List files
# star_files <- list.files("./results/STAR", 
#                          pattern="Aligned.sortedByCoord.out.bam$", 
#                          recursive = TRUE,
#                          full.names = TRUE)
# star_files <- c(snakemake@input$P635L_alignment_files, 
#                 snakemake@input$R636Q_alignment_files)
star_files <- unlist(snakemake@input)

# Create counts from STAR alignment
# featureCounts() from the Rsubread package
featurecounts_allsamples <- featureCounts(file = star_files,
                                          annot.ext = "resources/gencode.vM29.annotation.gtf",
                                          isGTFAnnotationFile = TRUE,
                                          strandSpecific = snakemake@params$strandness,
                                          isPairedEnd = TRUE,
                                          nthreads = snakemake@threads)
# Save object
saveRDS(featurecounts_allsamples, file = unlist(snakemake@output))