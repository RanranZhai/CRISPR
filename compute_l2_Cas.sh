DIR="/home/ranran/ldsc/1000G_EUR_Phase3_plink/"

cd $DIR

files=$(ls *.bim)

for i in $files;do
	var=${i##1000G.EUR.QC.}
	var=${var%%.bim}
	bim=/home/ranran/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.$var
	annot=/opt/working/projects/prj_026_CRISPR/Rev2203/Cas/$1/annotation/$var.annot.gz ## prefix of annot.gz file, without chromosome number, which was refer to $var
	out=/opt/working/projects/prj_026_CRISPR/Rev2203/Cas/$1/annotation/$var
	snp=/home/ranran/ldsc/listHM3.txt
	echo $var
	echo $bim
	echo $annot
	python /home/ranran/ldsc/ldsc.py --l2 --bfile $bim --ld-wind-cm 1 --annot $annot --out $out --print-snps $snp --thin-annot
done






