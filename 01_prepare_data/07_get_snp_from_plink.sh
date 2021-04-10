
#!/bin/bash

dir="/Users/anastasiia_hry/bio/datasets/snps"
dir="/home/ahryhorzhevska/mpip/datasets/2020_DexStim_Array_Human/snps"

# recode the Dex_genoData into .gen format + .sample
$dir/plink --bfile $dir/Dex_genoData_SNPs --recode oxford --out $dir/Dex_genoData_SNPs
echo "recoded Dex_genoData"

# retrieve the SNP position - consisting of chromosome, ID and position
cut -f 2 -d' ' $dir/Dex_genoData_SNPs.gen | sed 's/ /_/g' > $dir/SNPposition.txt
echo "done with retrieving info about SNP postion"

# retrieve the info with three numbers describing one sample (1st: homozygote reference, 2nd: heterozygote, 3rd: homozygote alternative)
cut -f 6- -d' ' $dir/Dex_genoData_SNPs.gen > $dir/plinkgen_X.txt
echo "done with retrieving info about samples"