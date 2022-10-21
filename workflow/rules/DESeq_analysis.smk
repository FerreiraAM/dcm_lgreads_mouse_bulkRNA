# DESeq analysis


# Create count matrices for DESeq analysis with start alignment
rule count_matrices:
  input:
    lambda wildcards: expand("results/STAR/{{experiment}}/{sample}/Aligned.sortedByCoord.out.bam", sample = config["samples"][wildcards.experiment])
  conda:
    "envs/process_dcm_mouse_bulkRNA_Ranalysis.yml"
  params:
    strandness = lambda wildcards: config["strand"][wildcards.experiment]
  output:
    "results/robject_featurecounts_{experiment}.rds"
  resources:
    mem_mb = "10G",
    time = "2-00:00:00"
  threads:
    6
  script:
    "scripts/DESeq_analysis_withSTARalignment_createcountmatrices.R"
    
