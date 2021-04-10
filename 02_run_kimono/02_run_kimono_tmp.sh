#!/bin/bash

src_dir=/home/ahryhorzhevska/mpip/code/dex-stim-human-kimono

# module load R

partition=pe
memory=200G
job_name=kimono

sbatch --job-name=$job_name --part=$partition --mem=$memory \
--wrap="Rscript --vanilla $src_dir/02_run_kimono/02_run_kimono_tmp.R"
