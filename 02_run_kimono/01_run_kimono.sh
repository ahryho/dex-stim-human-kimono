#!/bin/bash

src_dir=/home/ahryhorzhevska/mpip/code/dex-stim-human-kimono

# module load R

partition=pe
node=5
memory=900G
job_name=kimono_

treatment="dex"
sbatch --job-name=$job_name$treatment --part=$partition --nodelist=$partition$node --mem=$memory \
	--wrap="Rscript --vanilla $src_dir/02_run_kimono/01_run_kimono.R $treatment"

node=7
memory=900G
treatment="veh"
sbatch --job-name=$job_name$treatment --part=$partition --nodelist=$partition$node --mem=$memory \
	--wrap="Rscript --vanilla $src_dir/02_run_kimono/01_run_kimono.R $treatment"
