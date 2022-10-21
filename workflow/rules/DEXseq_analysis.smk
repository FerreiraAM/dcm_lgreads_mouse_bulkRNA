# Run DEXSeq

# Prepare the annotation
rule DEXSeq_annotation_preparation:
  input:
    gtf_file = "resources/gencode.vM29.annotation.gtf"
  output:
    "results/gencode_GRCm39/gencode.vM29.annotation.DEXSeq.chr.gff"
  conda:
    "envs/process_dcm_mouse_bulkRNA_Ranalysis.yml"
  shell:
    "python /g/steinmetz/ferreira/R-lib/4.1.2-foss-2021b/DEXSeq/python_scripts/dexseq_prepare_annotation.py {input.gtf_file} {output}"
    
# Sort bam files
rule sort_STAR_bam_files:
  input:
    alignment_file = "results/STAR/{mutation}/{sample}/Aligned.sortedByCoord.out.bam"
  output:
    "results/STAR/{mutation}/{sample}/Aligned.sortedByCoord.out_sorted.bam"
  conda:
    "envs/process_dcm_mouse_bulkRNA_MAPPING.yml"
  shell:
    "samtools sort -n {input.alignment_file} -o {output}"
    
# Count reads
rule DEXSeq_count_reads:
  input:
    gff_file = "results/gencode_GRCm39/gencode.vM29.annotation.DEXSeq.chr.gff",
    alignment_file = "results/STAR/{mutation}/{sample}/Aligned.sortedByCoord.out_sorted.bam"
  output:
    "results/STAR/{mutation}/{sample}/read_counts_sorted_reverse.txt"
  conda:
    "envs/process_dcm_mouse_bulkRNA_Ranalysis.yml"
  resources:
    time = "2-00:00:00"
  shell:
    "python /g/steinmetz/ferreira/R-lib/4.1.2-foss-2021b/DEXSeq/python_scripts/dexseq_count.py -p yes -f bam -s reverse {input.gff_file} {input.alignment_file} {output}"

# Clean reads_counts.txt files
# https://support.bioconductor.org/p/9143537/
rule clean_reads_counts:
  input:
    "results/STAR/{mutation}/{sample}/read_counts_sorted_reverse.txt"
  output:
    "results/STAR/{mutation}/{sample}/read_counts_sorted_reverse_clean.txt"
  conda:
    "envs/process_dcm_mouse_bulkRNA_MAPPING.yml"
  resources:
    time = "2-00:00:00"
  shell:
    "sed 's/\"//g' {input} > {output}"
    
# DEXSeq estimate dispersion
rule DEXSeq_dispersion:
  input:
    meta = config["samples_metadata"],
    P635L_files = expand("results/STAR/P635L/{sample}/read_counts_sorted_reverse_clean.txt", sample = config["samples"]["P635L"]),
    R636Q_files = expand("results/STAR/R636Q/{sample}/read_counts_sorted_reverse_clean.txt", sample = config["samples"]["R636Q"])
  conda:
    "envs/process_dcm_mouse_bulkRNA_Ranalysis.yml"
  output:
    "results/DEXSeq/robject_{mutation}_dexseq_estdisp_3conditions.rda",
    "results/DEXSeq/plot_{mutation}_dxd_dispersion_estimate_3conditions.png",
    "results/DEXSeq/robject_{mutation}_dexseq_estdisp_WTvsHOM.rda",
    "results/DEXSeq/plot_{mutation}_dxd_dispersion_estimate_WTvsHOM.png"
  resources:
    mem_mb = "40G",
    time = "2-00:00:00"
  threads:
    6
  script:
    "scripts/DEXSeq_analysis_dispersion.R"
    
# DEXSeq DEU test
rule DEXSeq_DEUtest:
  input:
    meta = config["samples_metadata"],
    files = expand("results/DEXSeq/robject_{mutation}_dexseq_estdisp_WTvsHOM.rda", mutation = config["mutation"])
  conda:
    "envs/process_dcm_mouse_bulkRNA_Ranalysis.yml"
  output:
    "results/DEXSeq/robject_{mutation}_dexseq_DEUres_WTvsHOM.rda"
  resources:
    mem_mb = "10G",
    time = "2-00:00:00"
  script:
    "scripts/DEXSeq_analysis_DEUtest.R"
