# Run STAR aligment

# Package
import glob

# infer input fastq files from dirs in config file and sample wildcard
def get_fastq_files(wildcards):
	indir = config["samples"][wildcards.mutation][wildcards.sample]
	file1 = glob.glob(indir + "/*" + wildcards.sample  + "_1_sequence.txt.gz")
	file2 = glob.glob(indir + "/*" + wildcards.sample  + "_2_sequence.txt.gz")
	return {"fastq1" : file1, "fastq2" : file2}

# Download gencode annotation
rule download_mouse_gencode_annotation:
  output: 
    annot = "resources/gencode.vM29.annotation.gtf"
  conda: 
    "../envs/process_dcm_mouse_bulkRNA_MAPPING.yml"
  params:
    url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M29/gencode.vM29.annotation.gtf.gz"
  shell:
    "wget -O {output.annot}.gz {params.url}; gunzip {output.annot}.gz"

# Donwload gencode fasta file
rule download_mouse_gencode_genome:
  output:
    genome = "resources/GRCm39.primary_assembly.genome.fa"
  conda:
    "../envs/process_dcm_mouse_bulkRNA_MAPPING.yml"
  params:
    url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M29/GRCm39.primary_assembly.genome.fa.gz"
  shell:
    "wget -O {output.genome}.gz {params.url}; gunzip {output.genome}.gz"

# Create STAR index
rule STAR_index:
  input:
    genome = "resources/GRCm39.primary_assembly.genome.fa",
    annot = "resources/gencode.vM29.annotation.gtf"
  output:
    dir = directory("resources/gencode_GRCm39/")
  conda:
    "../envs/process_dcm_mouse_bulkRNA_MAPPING.yml"
  resources:
    mem_mb = "40G"
  threads:
    6
  shell:
    "STAR --runThreadN 6 --runMode genomeGenerate --genomeDir {output.dir} --genomeFastaFiles {input.genome} --sjdbGTFfile {input.annot} --sjdbOverhang 100"

# Align reads with STAR
rule STAR_mapping:
  input:
    unpack(get_fastq_files),
    dir = "resources/gencode_GRCm39/"
  params:
    prefix = "results/STAR/{mutation}/{sample}/"
  output:
    "results/STAR/{mutation}/{sample}/Aligned.sortedByCoord.out.bam"
  conda:
    "../envs/process_dcm_mouse_bulkRNA_MAPPING.yml"
  resources:
    mem_mb = "30G",
    time = "2-00:00:00"
  threads:
    6
  shell:
    "STAR --genomeDir {input.dir} --runThreadN 6 --readFilesIn {input.fastq1} {input.fastq2} --outFileNamePrefix {params.prefix} --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard  --readFilesCommand zcat"

# STAR report info
rule STAR_report_info:
  input:
    expand("results/STAR/P635L/{sample}/Log.final.out", sample = config["samples"]["P635L"]),
    expand("results/STAR/R636Q/{sample}/Log.final.out", sample = config["samples"]["R636Q"]),
    expand("results/STAR/P635L_after_base_editing/{sample}/Log.final.out", sample = config["samples"]["P635L_after_base_editing"]),
    expand("results/STAR/R636Q_after_base_editing/{sample}/Log.final.out", sample = config["samples"]["R636Q_after_base_editing"])
  conda:
    "../envs/process_dcm_mouse_bulkRNA_MAPPING.yml"
  output:
    "results/table_STARalignment_total_number_reads.txt"
  script:
    "scripts/align_report.R"    
