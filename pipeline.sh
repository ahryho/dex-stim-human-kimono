#!/bin/bash
#
#SBATCH --job-name=kimono
#SBATCH --output=err_out/out_%A_%a.out
#SBATCH --error=err_out/err_%A_%a.err
#SBATCH --array=0-1
#SBATCH --mem=2G
#SBATCH --part=pe

src_data_methyl_dir_pre='/Users/anastasiia_hry/bio/datasets/methylation/20_DMA/'
data_kimono_dir='~/bio/datasets/kimono/'
wd_methyl_pre='/Users/anastasiia_hry/bio/code/mpip/dex-stim-human-differential-methyl-analysis'
wd_kimono='/Users/anastasiia_hry/bio/code/mpip/dex-stim-human-kimono/01_prepare_data'

cd $wd_methyl_pre

data=$src_data_methyl_dir_pre/"02_dmp/dmps_pval_adj_bcc_pcs.txt"
model='bcc_pcs'
delta_beta=(0.1)
pval_thrsh=(0.05)

for ibeta in ${delta_beta[@]}
do
    for ipval in ${pval_thrsh[@]}
    do
        dmp_sign=(("dmps_significant_with_mdl_"$model"_beta_"(expr $ibeta * 100)"_p_"(expr $ipval * 100)".txt"))
        section="get significant DMPs from the model "$model" with delta beta = "$ibeta" and p-value = "$ipval
        echo $section
        Rscript --vanilla $wd_methyl_pre/02_dmp/02_get_sign_dmps_skm.R --directory=$src_data_methyl_dir_pre --data=$data --model=$model

        section="prepare methylation data and prior 'CpG <-> Gene' for kimono"
        echo $section
        Rscript --vanilla $wd_kimono/ 03_b_prior_beta_mtrx.R  --directory=$src_data_methyl_dir_pre \
                                                              --data=\
                                                              --betamtrx="/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/10_final_qc_data/dex_methyl_beta_combat_mtrx.rds" \
                                                              --mapping="/mapping/mapping_cpg_gene_ensg.csv" \
                                                              --out=$data_kimono_dir

    done
done




# #BSUB -J combp[4-25]
# #BSUB -R rusage[mem=6]
# #BSUB -e logs/combp.%I.%J.err
# #BSUB -o logs/combp.%I.%J.out
# #BSUB -n 12
# #BSUB -R span[hosts=1]

# mkdir -p logs/ 
# mkdir -p data
# set -e

# #H=/home/brentp/src/combined-pvalues/examples/charm
# #cd $H

# SECTION=$1
# SECTION=COMB
# echo $SECTION

# if [ "$SECTION" = "GET" ]; then
# cd data
# # remove extra space and save.
# wget -O - http://rafalab.jhsph.edu/data/shores/natgen2009.csv \
#     | perl -pe 's/,\s/,/' > natgen2009.csv

# mkdir -p xys/ 
# cd xys/
# wget http://rafalab.jhsph.edu/data/shores/natgen2009.tgz
# tar xzvf natgen2009.tgz
# exit;
# fi

# if [ "$SECTION" = "NORMALIZE" ]; then

# R --slave --args data/natgen2009.csv data/methp.txt \
#     SampleID quantile < scripts/charm.normalize.R
# sed -i "1s/^/ID\t/" data/methp.txt
# exit;

# fi


# if [ "$SECTION" = "FIT" ]; then
# mkdir -p data/fit
# NCOL=2162407
# STEP=8000

# # get rid of errant space that messes up matching.
# for i in `awk -v cols=$NCOL -v step=$STEP 'BEGIN{for(i=2;i<cols;i+=step){print i }}'`; do
#     start=$i
#     end=$(($i + $STEP - 1))
#     nf=data/fit/$start-$end.split
#     bedp=data/fit/$start-$end.bed
#     #cut -f 1,$start-$end data/methp.txt > $nf
#     echo "~/local/bin/R --slave < scripts/fit.lm.R --args $nf > $bedp" \
#         | bsub  \
#         -e logs/$start-$end.err \
#         -o logs/$start-$end.out \
#         -R "rusage[mem=6]"
# done
# exit;
# fi


# if [ "$SECTION" = "MERGE" ]; then
# cat data/fit/*.bed | awk 'NR == 1 || $1 != "#chrom"' \
#     | sort -k1,1 -k2,2n > data/pvalues.bed
# exit;
# fi


# if [ "$SECTION" = "COMB" ]; then

# COLS=(fake chrom start end p.disease p.tissue)

# COL=${COLS[$LSB_JOBINDEX]}

# if [[ $LSB_JOBINDEX -gt 5 ]]; then
#     COL="p.disease-$LSB_JOBINDEX"
# fi

# if [[ $LSB_JOBINDEX -gt 25 ]]; then
#     COL="p.tissue-$LSB_JOBINDEX"
# fi

# PRE=data/quantile/$COL/$COL
# mkdir -p data/quantile/$COL

# <<DONE
# comb-p pipeline \
#     -c $LSB_JOBINDEX \
#     -s \
#     --seed 0.0005 \
#     --dist 80 --step 40 \
#     -p $PRE \
#     data/pvalues.bed

#     awk 'NR == 1 || ($5 > 5 && $7 < 0.001)' ${PRE}.regions-p.bed \
#             > ${PRE}.sig.regions.bed
# DONE
# # to check for # false positives on just the raw p-values

# # only print out the FDR column
# comb-p fdr -c $LSB_JOBINDEX data/pvalues.bed | awk 'BEGIN{OFS=FS="\t"}{ print $1,$2,$3,$NF }' > $PRE.fdr
# comb-p peaks -c $LSB_JOBINDEX --seed 0.0005 --dist 80 data/pvalues.bed > $PRE.peaks;
# comb-p peaks -c 4 --seed 0.0005 --dist 80 $PRE.fdr > $PRE.fdr.peaks;

# exit;

# fi

# if [ "$SECTION" = "SPLIT_FIT" ]; then

# mkdir -p data/split_fit
# mkdir -p logs/split_fit

# for start in $(seq 2 8000 2058592); do

#     f=data/split_fit/f$start.bed
#     end=$((start + 8000))
#     awk -vstart=$start -v end=$end 'NR == 1 || (NR > start && NR < end)' data/pvalues.bed | cut -f 1-4 > $f
#     echo "comb-p pipeline -c 4 -s --seed 0.0005 --dist 80 --step 40 -p data/split_fit/s$(basename $f .bed) $f" | bsub -J $(basename $f) \
#             -e logs/split_fit/$(basename $f .bed).err \
#             -o logs/split_fit/$(basename $f .bed).out \
#             -R "rusage[mem=4]"
#     echo $f
# done
# exit;

# fi

# if [ "$SECTION" = "PEAKS_RAW" ]; then
#     f=data/pvalues.bed
#     echo "
#     comb-p peaks -c 6 --seed 0.0005 --dist 80 $f > ${f/bed/peaks};
#     comb-p region_p -r ${f/bed/peaks} -p $f -s 40 -c 6 > data/raw.peaks.bed;
#     "| bsub -J comb-p-peaks -e logs/comb-peaks.err -o logs/comb-peaks.out
#     awk '$7 < 1' data/raw.peaks.bed # none !
#     exit;
# fi

# echo "send in a section to run, e.g.: bash run.sh FIT"
# echo "sections should be run in order:
#    GET
#    NORMALIZE
#    FIT
#    MERGE
#    COMB"



