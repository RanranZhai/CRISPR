
##CRISPR_effeciency ldsc

FIRST step: make ucsc .bed file from genomic region information

To run this script, 4 arguments will be needed
 
 Args[1] : file name of input, 
	   which needs three columns 
	   column 1 indicating chromsome (e.g. 1)
		 column 2 start position
		 column 3 end position (optional)

 Args[2] : number, indicating bed file start from original start
 
 Args[3] : number, indicating bed file end from original start

 Args[4] : file name of output (result)


## Rscirpt make_ucsc_bed.R Args[1] ... Args[4]
