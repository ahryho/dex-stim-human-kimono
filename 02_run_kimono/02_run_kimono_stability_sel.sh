#!/bin/bash
#
#SBATCH --job-name=kimono_diffGRN
#SBATCH --output=err_out/kimono_diffgrn_%A_%a.out
#SBATCH --error=err_out/kimono_diffgrn_%A_%a.err
#SBATCH --array=0-79
#SBATCH --mem=200Gb
#SBATCH --cpus-per-task=12
#SBATCH --partition=pe

export LC_CTYPE="en_US.UTF-8"
export LC_ALL="en_US.UTF-8"

treatment=("veh" "dex")
startnodes=($(seq 1 10000 400000))

nstartnodes=${#startnodes[@]}
ndex=${#treatment[@]}

#get region and dex index for each job id
istartnode=$((SLURM_ARRAY_TASK_ID / ndex)) #divide task id by number of dex status
idex=$(($SLURM_ARRAY_TASK_ID%$ndex))

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "${startnodes[$istartnode]}"
echo "${treatment[$idex]}"
echo $istartnode
echo $idex

Rscript --vanilla 02_run_kimono_stability_sel.R ${treatment[$idex]} ${startnodes[$istartnode]}