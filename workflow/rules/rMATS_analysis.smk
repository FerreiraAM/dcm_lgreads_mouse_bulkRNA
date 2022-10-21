# rMATS

# rule splicing_event_with_rMATS:
#   input:
#     gtf_file = "resources/gencode.vM29.annotation.gtf",
#     wt_file = "results/rmats_directory/P635L_WT_files.txt",
#     het_file = "results/rmats_directory/P635L_HET_files.txt"
#   output:
#     "A3SS.MATS.JC.txt",
#     "A5SS.MATS.JC.txt",
#     "MXE.MATS.JC.txt",
#     "RI.MATS.JC.txt",
#     "SE.MATS.JC.txt",
#     dir = "results/rmats_directory/P635L_HETvsWT/"
#   conda:
#     "envs/process_dcm_mouse_bulkRNA_rMATS.yml"
#   resources:
#     mem_mb = "10G",
#     time = "1-00:00:00"
#   params:
#     readLength = 160,
#     typeread = "paired"
#   threads:
#     6
#   shell:
#     "python /g/steinmetz/ferreira/software/miniconda3/envs/rMATS_env/rMATS/rmats.py --b1 {input.het_file} --b2 {input.wt_file} --gtf {input.gtf_file} -t {params.typeread} --readLength {params.readLength} --nthread {threads} --od {output.dir}"
# conda activate rMATS_env

# # rMATS data analysis - overlap of SJ
# rule rMATS_overlapSJ:
#   conda:
#     "envs/process_dcm_mouse_bulkRNA_Ranalysis.yml"
#   output:
#     "rmats_directory/ls_sameevent_pergene.rds"
#   resources:
#     mem_mb = "10G",
#     time = "2-00:00:00"
#   script:
#     "scripts/rmats_data_analysis_cluster.R"
