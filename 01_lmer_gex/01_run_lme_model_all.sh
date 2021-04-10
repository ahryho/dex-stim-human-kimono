#!/bin/bash

src_dir=/home/ahryhorzhevska/mpip/code/dex-stim-human-kimono/01_lmer_gex
rslt_dir=/binder/mgp/datasets/2020_DexStim_Array_Human/gene_expression/20_DEA/01_lme_models/

rslt_lme=$rslt_dir/lmer_all_plus_bcc.csv
out_fn=lmem_gex_all.out
# err_fn=.err

# module load R

partition=pe
memory=100G
job_name=lme_gex

sbatch --job-name=$job_name --part=$partition --mem=$memory --output=$out_fn --wrap="Rscript --vanilla $src_dir/01_lme_model_all.R $rslt_lme"
