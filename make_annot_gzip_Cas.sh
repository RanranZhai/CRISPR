
#!/bin/bash

DIR="/home/ranran/ldsc/1000G_EUR_Phase3_plink/"

cd $DIR

files=$(ls *.bim)

echo $files

mkdir /opt/working/projects/prj_026_CRISPR/Cas/$1/annotation

for i in $files;do
	var=${i##1000G.EUR.QC.}
	var=${var%%.bim}
	bim=/home/ranran/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.$var.bim
	annot=/opt/working/projects/prj_026_CRISPR/Cas/$1/annotation/$var.annot
	bed=/opt/working/projects/prj_026_CRISPR/Cas/$1/Top10.bed
	echo $var
	echo $bim
	echo $annot
	echo $bed
	python /home/ranran/ldsc/make_annot.py  --bed-file $bed --bimfile $bim --annot-file $annot 
done

### first argument : ucsc .bed file (e.g. )
### second argument : explaintory (e.g. TTTA, TTTG or TTTC)


gzip -r /opt/working/projects/prj_026_CRISPR/Cas/$1/annotation/*.annot



