#!/bin/bash

src_dir=/home/ahryhorzhevska/mpip/code/dex-stim-human-kimono

# module load R

partition=pe
node=7
memory=200G
job_name=kimono_veh

sbatch --job-name=$job_name --part=$partition --nodelist=$partition$node --mem=$memory \
        --wrap="Rscript --vanilla $src_dir/02_run_kimono/01_run_kimono_veh.R"
