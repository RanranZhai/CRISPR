
#!/bin/bash

DIR="ldsc/1000G_EUR_Phase3_plink/"

cd $DIR

files=$(ls *.bim)

echo $files

for i in $files;do
	var=${i##1000G.EUR.QC.}
	var=${var%%.bim}
	bim=/Users/wenzheng/projects/prj_ldsc_final/1000G_EUR_Phase3_plink/1000G.EUR.QC.$var.bim
	annot=/Users/ranran/CRISPR_annotation/TTTN/TTTA_4bp.$var.annot.gz
	bed=/Users/ranran/CRISPR_effeciency/TTTN/TTTA_4bp.$var.bed
	echo $var
	echo $bim
	echo $annot
	echo $bed
	python /Users/ranran/ldsc/make_annot.py  --bed-file $bed --bimfile $bim --annot-file $annot
done










