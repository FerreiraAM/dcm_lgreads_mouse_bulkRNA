# activate conda envrs
conda activate rMATS_env
# script
# python /g/steinmetz/ferreira/software/miniconda3/envs/rMATS_env/rMATS/rmats.py
# Command
# python /g/steinmetz/ferreira/software/miniconda3/envs/rMATS_env/rMATS/rmats.py --b1 /g/steinmetz/project/dcm_lgreads/mouse_bulkRNA/results/STAR/P635L/21s003516/Aligned.sortedByCoord.out_sorted.bam,/g/steinmetz/project/dcm_lgreads/mouse_bulkRNA/results/STAR/P635L/21s003517/Aligned.sortedByCoord.out_sorted.bam,/g/steinmetz/project/dcm_lgreads/mouse_bulkRNA/results/STAR/P635L/21s003518/Aligned.sortedByCoord.out_sorted.bam,/g/steinmetz/project/dcm_lgreads/mouse_bulkRNA/results/STAR/P635L/21s003519/Aligned.sortedByCoord.out_sorted.bam,/g/steinmetz/project/dcm_lgreads/mouse_bulkRNA/results/STAR/P635L/21s003520/Aligned.sortedByCoord.out_sorted.bam --b2 /g/steinmetz/project/dcm_lgreads/mouse_bulkRNA/results/STAR/P635L/21s003506/Aligned.sortedByCoord.out_sorted.bam,/g/steinmetz/project/dcm_lgreads/mouse_bulkRNA/results/STAR/P635L/21s003507/Aligned.sortedByCoord.out_sorted.bam,/g/steinmetz/project/dcm_lgreads/mouse_bulkRNA/results/STAR/P635L/21s003508/Aligned.sortedByCoord.out_sorted.bam,/g/steinmetz/project/dcm_lgreads/mouse_bulkRNA/results/STAR/P635L/21s003509/Aligned.sortedByCoord.out_sorted.bam,/g/steinmetz/project/dcm_lgreads/mouse_bulkRNA/results/STAR/P635L/21s003510/Aligned.sortedByCoord.out_sorted.bam --gtf resources/gencode.vM29.annotation.gtf -t paired --readLength 50 --nthread 1 --od rmats_directory --tmp rmats_directory/tmp_output
# P635L
# HOM vs WT
python /g/steinmetz/ferreira/software/miniconda3/envs/rMATS_env/rMATS/rmats.py --b1 /g/steinmetz/project/dcm_lgreads/mouse_bulkRNA/rmats_directory/P635L_HOM_files.txt --b2 /g/steinmetz/project/dcm_lgreads/mouse_bulkRNA/rmats_directory/P635L_WT_files.txt --gtf resources/gencode.vM29.annotation.gtf -t paired --readLength 160 --nthread 1 --od rmats_directory/P635L_HOMvsWT/
# HET vs WT
python /g/steinmetz/ferreira/software/miniconda3/envs/rMATS_env/rMATS/rmats.py --b1 /g/steinmetz/project/dcm_lgreads/mouse_bulkRNA/rmats_directory/P635L_HET_files.txt --b2 /g/steinmetz/project/dcm_lgreads/mouse_bulkRNA/rmats_directory/P635L_WT_files.txt --gtf resources/gencode.vM29.annotation.gtf -t paired --readLength 160 --nthread 1 --od rmats_directory/P635L_HETvsWT/
# R636Q
# HOM vs WT
python /g/steinmetz/ferreira/software/miniconda3/envs/rMATS_env/rMATS/rmats.py --b1 /g/steinmetz/project/dcm_lgreads/mouse_bulkRNA/rmats_directory/R636Q_HOM_files.txt --b2 /g/steinmetz/project/dcm_lgreads/mouse_bulkRNA/rmats_directory/R636Q_WT_files.txt --gtf resources/gencode.vM29.annotation.gtf -t paired --readLength 160 --nthread 1 --od rmats_directory/R636Q_HOMvsWT/
# HET vs WT
python /g/steinmetz/ferreira/software/miniconda3/envs/rMATS_env/rMATS/rmats.py --b1 /g/steinmetz/project/dcm_lgreads/mouse_bulkRNA/rmats_directory/R636Q_HET_files.txt --b2 /g/steinmetz/project/dcm_lgreads/mouse_bulkRNA/rmats_directory/R636Q_WT_files.txt --gtf resources/gencode.vM29.annotation.gtf -t paired --readLength 160 --nthread 1 --od rmats_directory/R636Q_HETvsWT/

# Summary file available for older version
#python /g/steinmetz/ferreira/software/miniconda3/envs/rMATS_env/rMATS_P/summary.py 
