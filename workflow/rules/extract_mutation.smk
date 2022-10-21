# Extract reads at the mutation's loci to visualize in IGV

# Create bam index
rule create_bam_idx:
  input:
    "results/STAR/{mutation}/{sample}/Aligned.sortedByCoord.out.bam"
  conda:
    "envs/process_dcm_mouse_bulkRNA_MAPPING.yml"
  output:
    "results/STAR/{mutation}/{sample}/Aligned.sortedByCoord.out.bam.bai"
  shell:
    "samtools index {input}"

# Extract regions of mutations  
rule extract_mutation_region:
  input:
     "results/STAR/{experiment}/{sample}/Aligned.sortedByCoord.out.bam",
     "results/STAR/{experiment}/{sample}/Aligned.sortedByCoord.out.bam.bai"
  conda:
    "envs/process_dcm_mouse_bulkRNA_MAPPING.yml"
  output:
    "results/STAR/IGV_mutation_check/{experiment}_{sample}_rbm20_mutation.bam",
    "results/STAR/IGV_mutation_check/{experiment}_{sample}_rbm20_mutation.bam.bai"
  params:
    region = lambda wildcards: config["mutation_regions"][wildcards.experiment]
  shell:
    """
    samtools view -b {input} {params.region} > {output}
    samtools index {output}
    """
