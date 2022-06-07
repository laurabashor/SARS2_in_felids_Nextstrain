#!/bin/bash
output="/home/lbashor/cat_SARS2_manuscript_June2021/iVar/iVar_11_12_21/output/variant_tables/"

cat InputFileList2.txt | while read -r line
do

ivar filtervariants -t 0.5 -p ${output}${line}.filtered ${line}_ivar.tsv ${line}_R_ivar.tsv

done
