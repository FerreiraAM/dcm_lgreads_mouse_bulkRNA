# ALL rules for kallisto alignment

# All rule alignment
rule all_aligment:
  input:
    expand("results/P635L/{sample}/abundance.tsv", sample = config["samples"]["P635L"]),
    expand("results/R636Q/{sample}/abundance.tsv", sample = config["samples"]["R636Q"])
    
# Download index
rule download_genome:
	output: 
		idx = "resources/tx_indices/Mus_musculus.GRCm38.cdna.all.release-94_k31.idx"
	conda: 
		"envs/process_dcm_mouse_bulkRNA_MAPPING.yml"
	params:
		url = "https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/94/Mus_musculus.GRCm38.cdna.all.release-94_k31.idx.gz"
	shell:
		"wget -O {output.idx}.gz {params.url}; gunzip {output.idx}.gz"
		
# Run mapping
rule run_kallisto:
	input: 
	  unpack(get_fastq_files),
		index = "resources/tx_indices/Mus_musculus.GRCm38.cdna.all.release-94_k31.idx"
	output:
		"results/{mutation}/{sample}/abundance.tsv",
		"results/{mutation}/{sample}/abundance.h5"
	params:
	  outdir = "results/{mutation}/{sample}"
	conda: 
		"envs/process_dcm_mouse_bulkRNA_MAPPING.yml"
	shell:
		"kallisto quant -i {input.index} -o {params.outdir} {input.fastq1} {input.fastq2}"
		
# DESeq analysis with kallisto alignment    
rule deseq_analysis:
  input:
    meta = config["samples_metadata"],
    P635L_files = expand("results/P635L/{sample}/abundance.tsv", sample = config["samples"]["P635L"]),
    R636Q_files = expand("results/R636Q/{sample}/abundance.tsv", sample = config["samples"]["R636Q"])
  output:
    "results/total_nb_reads.csv",
    "results/plot_pca_all_mice_genes.pdf",
    "results/plot_pca_P635L_mice_genes.pdf",
    "results/plot_pca_R636Q_mice_genes.pdf",
    "results/table_res_P635L_homvswt_sig.txt",
    "results/table_res_P635L_hetvswt_sig.txt",
    "results/table_res_R636Q_homvswt_sig.txt",
    "results/table_res_R636Q_hetvswt_sig.txt",
    "results/heatmap_R636Q.pdf",
    "results/heatmap_P635L.pdf"
  conda:
    "envs/process_dcm_mouse_bulkRNA_Ranalysis.yml"
  script:
    "scripts/DESeq_analysis.R"
