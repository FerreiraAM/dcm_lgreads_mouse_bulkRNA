# Studying RBM20 mutations in mice

## Experiment descriptions

2 experiments were performed by Markus Grosch. The sequencing was performed using 
Illumina (bulk RNA-seq) from mouse hearts.

Two main mutations are studied in the RBM20 gene: P635L located at chr19:53831671 
(C wt > T mut) and R636Q located at chr19:53831674/675 (GT wt > AG mut).

### Induce mutation - experiment 1

Mutations were generated by CRISPR Cas9. Samples information are listed in 
config/samples_metadata.txt. There are 5 samples for each group: 
wild type (WT); heterozygote (HET); homozygote (HOM).

### After base editing - experiment 2

Multiple base editors were used. Samples information are listed in 
config/samples_afterbaseediting_metadata.txt.

## Analysis

The analysis was done using *Snakemake* [Mölder et al. 2021] and followed the *Snakemake* 
recommended directory structure (except for the *rMATS* analysis [Shen et al. 2014]).

### Alignment

#### STAR

The alignment of the different samples is performed using *STAR* [Dobin et al. 2013]. 
We use the GENCODE [Frankish et al. 2021] mouse annotation version vM29 with 
the primary assembly GRCm39 genome. 
We create the indexes and then align the reads for each sample using the default options
of the *STAR* aligner.

Script located at:

- workflow/rules/STAR_alignment.smk

### IVG visualization

We extract the sequences from the bam file to visualize the mutation's loci in 
*IGV* [Robinson et al. 2017].

Script located at:

- workflow/rules/extract_mutation.smk

### DESeq2

We use DESeq2 [Love et al. 2019] to perform the differential expression analysis 
using the R programming language [R Core Team (2022)]. We perform each comparison 
(listed below) for each mutation/conditions associated with each experiment using the count 
matrices that were created from the BAM files thanks to the Rsubread R package 
[Liao et al. 2019]. We created the DESeq2 object and use the main function (DESeq) 
that performs a default analysis through the following steps described in 
[Love et al. 2019]: (i) estimation of the size factors; (ii) estimation of the 
dispersion; (iii) negative binomial Generalized Linear Model fitting and Wald statistics.

Comparison list:

- P635L: HET vs WT; HOM vs WT;
- R636Q: HET vs WT; HOM vs WT;
- P635L after base editing: Nterm_SpRY_and_Cterm_SpRY_gRNA5 vs Nterm_NRTH_Abe8e_and_Cterm_gRNA5;
Nterm_SpRY_and_Cterm_SpRY_gRNA5 vs PBS; Nterm_SpRY_and_Cterm_SpRY_gRNA5 vs none; 
Nterm_NRTH_Abe8e_and_Cterm_gRNA5 vs PBS; Nterm_NRTH_Abe8e_and_Cterm_gRNA5 vs none; 
PBS vs none;
- R636Q after base editing: Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6 vs PBS; 
Nterm_NRTH_Abe8e_and_Cterm_NRCH_gRNA6 vs none; PBS vs none.

We computed the log2 fold change (log2FC) per individual for each mutation, 
using the average per gene from the WT samples of one experiment’s mutation: 
*log2FC for gene A = log2(value of gene A / WT average for gene A)*.

For visualization, we plotted different heatmaps with the top genes expression 
as well as with the log2FC values.

Scripts located at:

- workflow/rules/DESeq_analysis.smk: to create the count matrices;
- workflow/scripts/DESeq_analysis_withSTARalignment.R: to perform the DESeq2 
analysis and plot all heatmaps.

### rMATS

We use *rMATS* [Shen et al. 2014] to detect differential alternative splicing events. 
Unfortunately, *rMATS* was not easily compatible with snakemake at the time of the analysis. 
We then created its own conda environment that can be activated using *conda activate rMATS_env*.

Using the python script provided by the software, we performed the different 
comparisons for the samples that were base edited. We used all the samples of 
one condition for one mutation as replicated and specified that the sequences 
were paired and the read length was 160 base pairs (information extracted from the 
STAR reports). We imported the results in R and analyzed the results from the 
Junctions Counts (JC) files. We identified splice junction events that were 
overlapping between different conditions and filtered for significant events. 
An event was considered significant if the False Discovery Rate (FDR) was inferior 
to 0.01 and the Percentage Spliced in (PSI) value was either superior to 0.1 or 
inferior to -0.1. We classified these significant events in three categories: 
rescued, mis-spliced or unchanged. We computed the average PSI difference between 
untreated samples and the base-edited samples, and considered the absolute 
difference. We classified the events as unchanged if 
the absolute difference of PSI was inferior than 0.1. The remaining events were 
either classified as rescued or mis-spliced. In addition, we defined the PSI 
values of untreated samples as the original value x_original, the PSI values 
of base-edited samples as the edited value x_edited and used the following criteria:

- Rescued:
x_original > 0 and x_edited >= 0 or -0.2 <= x_edited <= 0.2 and x_original > x_edited;
x_original < 0 and x_edited <= 0 or -0.2 <= x_edited <= 0.2 and x_original < x_edited

- Mis-spliced:
x_original > 0 and x_edited >= 0 and x_original < x_edited;
x_original < 0 and x_edited >= 0 or 0.2 <= x_edited >= 0.2 and x_original > x_edited;
x_original > 0 and x_edited -0.2;
x_original < 0 and x_edited > 0.2

All *rMATS* results are located outside the *results* directory, in its own directory:
*rmats_directory/*.

Scripts located at:

- rmats_directory/README.txt: list of the commands to run the *rMATS* comparisons. The input files used to perform the comparison were created manually (*mutation_condition/genotype_files.txt*);
- workflow/scripts/rmat_data_analysis.R: R scripts analyzing the output of *rMATS* python script.

## Appendix

### Kallisto Alignment

The initial alignment was performed using the *kallisto* mapper [Bray et al. 2016]. 
However, we then switched to STAR as it was a requirement for DEXseq.

Script located at:

- rules/kallisto_alignment.smk

### DEXseq

The initial splice junction events analysis was performed using *DEXseq* R package 
[Anders et al. 2012].

Script located at:

- rules/DEXSeq_analysis.smk
