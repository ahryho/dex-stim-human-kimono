#!/bin/bash

src_dir=/home/ahryhorzhevska/mpip/code/dex-stim-human-kimono
data_dir=/binder/mgp/datasets/2020_DexStim_Array_Human/eQTM/
treatment="veh"

# module load R

partition=pe
node=4
memory=200G
job_name=eqtm_$treatment

sbatch --job-name=$job_name --part=$partition --nodelist=$partition$node --mem=$memory \
	--wrap="Rscript --vanilla $src_dir/04_eqtm/03_run_eqtm.R $treatment $data_dir"