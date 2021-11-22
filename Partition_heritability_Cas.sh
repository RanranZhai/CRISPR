## run Cas annotation using baselineV2.2
##

DIR='/home/ranran/ldsc/formatted_gwas_29traits/'

cd $DIR

files=$(ls *.sumstats.gz)

#mkdir /Users/ranran/Projects/CRISPR_new/$1/effeciency
mkdir /opt/working/projects/prj_026_CRISPR/Cas/$1/results
for file in $files;do
	annot=/opt/working/projects/prj_026_CRISPR/Cas/$1/annotation/
	base=/home/ranran/ldsc/baseline_v2.2/baselineLD.
	out=/opt/working/projects/prj_026_CRISPR/Cas/$1/results/$file
	weight=/home/ranran/ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
	frq=/home/ranran/ldsc/1000G_Phase3_frq/1000G.EUR.QC.
	echo $file
	python /home/ranran/ldsc/ldsc.py --h2 $file --ref-ld-chr $annot,$base --w-ld-chr $weight --overlap-annot --frqfile-chr $frq --out $out
done


