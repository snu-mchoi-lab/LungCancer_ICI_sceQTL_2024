#!/bin/bash

out_prefix="$1"
path_plink="/home/rstudio/miniconda3/bin/plink"
path_original="/path/to/geno/230628_chrALL_sample_ALLIO" 
temp_fname="/home/rstudio/.local/temp/temp_plink"
out_fname="plink_subset/${out_prefix}"

$path_plink --bfile $path_original \
  --keep-allele-order \
  --keep ${temp_fname}.subset_ids.txt \
  --update-ids ${temp_fname}.update_ids.txt \
  --make-bed --out ${temp_fname}.update_ids
  
$path_plink --bfile ${temp_fname}.update_ids \
  --keep-allele-order \
  --update-sex ${temp_fname}.update_sex.txt \
  --make-bed --out ${temp_fname}.update_sex
  
$path_plink --bfile ${temp_fname}.update_sex \
  --maf 0.01 --hwe 1e-5 --geno 0.03 --mind 0.03 \
  --output-chr chrM \
  --keep-allele-order \
  --make-bed --pca --out ${out_fname}_maf01_HWEe5_geno03_mind03 && \
  rm ${temp_fname}*
