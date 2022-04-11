## OR for Cas annotation and functiona annotations

flb <- list.files('/home/ranran/ldsc/baseline_v2.2')
flb <- flb[grep('annot.gz', flb)]

baseline <- c()
for (i in 1:length(flb)) {
	tmp <- read.table(gzfile(paste0('/home/ranran/ldsc/baseline_v2.2/', flb[i])), header = T, stringsAsFactors = F)
	baseline <- rbind(baseline, tmp)
	cat(i, ' ')
}


## get interger baseline 
cla <- c()
for (i in 6:ncol(baseline)) {
	tmp <- class(baseline[,i])
	cla <- c(cla, tmp)
}

idx <- which(cla == 'integer')

base <- baseline[,-c(1:5)]
base <- base[,idx]
base <- base[,-82]

## remove flanking annotation
idx <- grep('flanking.500', colnames(base))
base <- base[,-idx]

## remove MAFbin annotation
idx <- grep('MAFbin', colnames(base))
base <- base[,-idx]




dir <- list.dirs('/opt/working/projects/prj_026_CRISPR/Rev2203/Cas/')
#dir <- dir[grep('top10', dir, fixed = T)]
dir <- dir[grep('annotation', dir, fixed = T)]


pams <- unlist(strsplit(dir, '/opt/working/projects/prj_026_CRISPR/Rev2203/Cas//'))[seq(2, 2*length(dir),2)]
pams <- unlist(strsplit(pams, '/annotation'))
#pams[length(pams)] <- 'GC_content'


OR <- matrix(NA, ncol = ncol(base), nrow = length(dir))
#P <- matrix(NA, ncol = ncol(base), nrow = length(dir))

for (i in 1:length(dir)) {
	fl <- list.files(dir[i])
	fl <- fl[grep('annot.gz', fl)]
	dd <- c()
	for (j in 1:length(fl)) {
		tmp <- read.table(paste0(dir[i], '/', fl[j]), header = T, stringsAsFactors = F)
		dd <- rbind(dd, tmp)
		cat(j, ' ')
	}
	for (j in 1:ncol(base)) {
		OR[i,j] <- fisher.test(table(base[, j], dd[,1]))$estimate
		#P[i,j] <- fisher.test(table(base[, j], dd[,1]))$p
		cat(j, ' ')
	}
	cat(i, '\n')
}

colnames(OR) <- colnames(base)
rownames(OR) <- pams

saveRDS(OR, file = '/opt/working/projects/prj_026_CRISPR/Rev2203/Cas_SNP_baseline_OR.rds')



OR <- readRDS('~/Documents/Projects/NGG/Revision2203-debug/Cas_SNP_baseline_OR.rds')

new_names <- c('Coding (UCSC)' ,'Conserved (LindbladToh)' ,'CTCF (Hoffman)' ,'DGF (ENCODE)' ,'DHS Peak (Trynka)' ,'DHS (Trynka)'
	,'Enhancer (Andersson)' ,'Enhancer (Hoffman)' ,'Fetal DHS (Trynka)' ,'H3K27ac (Hnisz)' ,'H3K27ac (PGC2)' ,'H3K4me1 Peaks (Trynka)' 
	,'H3K4me1 (Trynka)' ,'H3K4me3 Peaks (Trynka)' ,'H3K4me3 (Trynka)' ,'H3K9ac Peaks (Trynka)' ,'H3K9ac (Trynka)' ,'Intron (UCSC)'
	,'Promoter Flanking (Hoffman)' ,'Promoter (UCSC)' ,'Repressed (Hoffman)' ,'H3K27ac (Hnisz)' ,'TFBS (ENCODE)' ,'Transcribed (Hoffman)' ,'TSS (Hoffman)' 
	,'3\' UTR (UCSC)' ,'5\' UTR (UCSC)' ,'Weak Enhancer (Hoffman)' ,'GERP RS > 4' ,'synonymous' ,'non_synonymous' ,'Conserved (Vertebrate)'
	,'Conserved (Mammal)' ,'Conserved (Primate)' ,'BivFlnk' ,'Promoter (Villar)' ,'Enhancer (Villar)' ,'Ancient Promoter' ,'Ancient Enhancer' ,'Promoter of ExAC genes')

colnames(OR) <- new_names


## determine Cas order
cas_info <- read.csv('~/Documents/Projects/NGG/Revision-1/CRISPR_PAM_full_table2201.csv', header = T, stringsAsFactors = F)
cas <- unique(cas_info$Cas)


or <- OR[cas,]

colMain <- colorRamp2(c(-4, 0, 2), c("#BC3C29", "white", "#0072B5"))

M = 0.7
pdf('~/Documents/Projects/NGG/Revision2203-debug/Cas_annotation_OR_2203.pdf',16,8)
Heatmap(as.matrix(log(or)), name = 'Enrichment', col = colMain, 
		show_heatmap_legend = T,
		heatmap_legend_param = list(
				title = "log(OR)", at = c(-4, -2, 0, 2)#, 
				#labels = c("neg_two", "zero", "pos_two")
				),
		height = unit(21*M, 'cm'), width = unit(40*M, 'cm'),
		#cell_fun = Pval_fun,
		row_order = 1:21, column_order = 1:40, 
		#row_split = row_splits, 
		#column_split = col_splits, 
		#column_title = NULL, 
		#row_title_gp = gpar(fontsize = 13, col = trait_cate_col), row_title_rot = 90, row_title_side = 'right',
		#bottom_annotation = featr_anoot,
		row_names_gp = gpar(fontsize = 12), #col = cas_col), 
		row_names_side = 'left', column_names_side = 'top',
		column_names_gp = gpar(fontsize = 12), column_names_rot = 60
		)
dev.off()













