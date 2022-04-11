### This script is to count all CRISPR-Cas PAM combination (from frontiers in genetics table) paper in 
## https://www.frontiersin.org/articles/10.3389/fgene.2020.00851/full
## ----------------------------the CRISPR full table was updated @2022.01.20-------------


#cas_info <- read.csv('~/Documents/Projects/NGG/Revision-1/CRISPR_PAM_full_table2201.csv', header = T, stringsAsFactors = F)
cas_info <- read.csv('/opt/working/projects/prj_026_CRISPR/CRISPR_PAM_full_table2201.csv', header = T, stringsAsFactors = F)
ddd <- cas_info[!duplicated(cas_info$PAM),]

pams <- ddd$PAM
r_pams <- ddd$rPAM
pam_len <- ddd$pam_len


idx <- which(nchar(pams) > 1)
pams <- pams[idx]
r_pams <- r_pams[idx]


for (i in 1:length(pams)) {
	#if (!file.exists(paste0('/opt/working/projects/prj_026_CRISPR/Rev2203/PAM_auto_2203/', pams[i], '_PAM_count_both_strand.txt'))) {
		#system(paste0('python /opt/working/projects/prj_026_CRISPR/fetch.py /mnt2/CRISPR/human_g1k_v37.fasta ', pams[i], ' > /mnt2/CRISPR/', pams[i], '_hg19_pos.tsv'))
		#system(paste0('python /opt/working/projects/prj_026_CRISPR/fetch.py /mnt2/CRISPR/human_g1k_v37.fasta ', r_pams[i], ' > /mnt2/CRISPR/', r_pams[i], '_hg19_pos.tsv'))
		# above two lines have been ran 
		#system(paste0('Rscript /opt/working/projects/prj_026_CRISPR/Rev2203/make_bed_PAM2203.R ', pams[i], ' ', r_pams[i]))
		system(paste0('Rscript /opt/working/projects/prj_026_CRISPR/Rev2203/make_bed_PAM0309.R ', pams[i], ' ', r_pams[i]))
	#}
	cat(pams[i], '\n')
}


### count GC content seperately (NG PAM also aplicable)
----- see script @/Users/lanez/Documents/Projects/NGG/Revision-1/Segments_non_N_length.R -----
### count GC content seperately (NG PAM also aplicable)





## annot cas protein
cas_info <- read.csv('/opt/working/projects/prj_026_CRISPR/CRISPR_PAM_full_table2201.csv', header = T, stringsAsFactors = F)
cas <- unique(cas_info$Cas)


## get chromosome length
load('/opt/ShaoYa/Users/ranran/PAM/GRC37_genome/CHRs_length.RData')
N <- 20000
M <- dd[1,23]/N
M <- round(M)


cut_length <- readRDS('/opt/working/projects/prj_026_CRISPR/human_g1k_v37_20000_cut_base_count.rds')
seg_length <- cut_length$non_N_length ## get segments non N lengths


for (Cas in cas) {
	dir.create(paste0('/opt/working/projects/prj_026_CRISPR/Rev2203/Cas/', Cas))
	PAM <- cas_info[which(cas_info$Cas == Cas), 'PAM']
	aa <- read.table(paste0('/opt/working/projects/prj_026_CRISPR/PAM/', PAM[1], '_PAM_count_both_strand.txt'), header = T, stringsAsFactors = F)
	aa <- aa[,1:3]
	for (pam in PAM) {
		tmp <- read.table(paste0('/opt/working/projects/prj_026_CRISPR/Rev2203/PAM_auto_2203/', pam, '_PAM_count_both_strand.txt'), header = T, stringsAsFactors = F)
		aa <- cbind(aa, tmp[,4])
	}
	colnames(aa)[-c(1:3)] <- PAM

	if (ncol(aa) > 4) {
		aa$PAM_sum <- rowSums(aa[,-c(1:3)])
	} else {
		aa$PAM_sum <- aa[,4]
	}
	#aa$PAM_sum2 <- aa$PAM_sum/seg_length

	yy <- aa[order(aa[,'PAM_sum'], decreasing = T),]
	top10 <- yy[1:2000,] ## get Top10%

	
	options(scipen = 200)
	Top10 <- top10[,c(1,2,3)]
	Top10$Chr <- paste0('chr', Top10$Chr)

	write.table(Top10, file = paste0('/opt/working/projects/prj_026_CRISPR/Rev2203/Cas/', Cas, '/Top10.bed'), col.names = F, row.names = F, quote = F, sep = '\t')
	system(paste0('bash /opt/working/projects/prj_026_CRISPR/Rev2203/make_annot_gzip_Cas.sh ', Cas))
	cat(Cas, 'done\n')
}

#write.table(Top10, file = paste0('/opt/working/projects/prj_026_CRISPR/Cas2/', 'AT_content', '/Top10.bed'), col.names = F, row.names = F, quote = F, sep = '\t')

## compute ld score
for (Cas in cas) {
	system(paste0('bash /opt/working/projects/prj_026_CRISPR/Rev2203/compute_l2_Cas.sh ', Cas))
}
## partition heritability
for (Cas in cas) {
	system(paste0('bash /opt/working/projects/prj_026_CRISPR/Rev2203/Partition_heritability_Cas.sh ', Cas))
}




## results


cas_info <- read.csv('/opt/working/projects/prj_026_CRISPR/CRISPR_PAM_full_table2201.csv', header = T, stringsAsFactors = F)
cas <- unique(cas_info$Cas)

names <- c('Alzheimer\'s disease', 'Body mass index', 'Breast cancer', 'Birth weight', 'Coronary artery disease', 'Bipolar disorder', 'Educational attainment', 'Celiac disease', 'Inflammatory bowel disease', 'Ulcerative colitis', 
	'Fasting glucose', 'Fasting insulin', 'Fasting proinsulin', 'HDL cholesterol', 'Height', 'LDL cholesterol', 'Lifespan', 'Type 2 diabetes', 'Major depression disorder', 'Attention deficit hyperactivity disorder', 'Autism spectrum disorder', 
	'Rheumatoid arthritis', 'Schizophrenia', 'Ever smoke', 'Total cholesterol', 'Triglycerides', 'Telomere length', 'Waist circumference', 'Waist-hip ratio')



res <- c()
for (Cas in cas) {
	dir <- paste0('/opt/working/projects/prj_026_CRISPR/Rev2203/Cas/', Cas, '/results/')
	fl <- list.files(dir)
	fl <- fl[grep('gz.results',fl)]

	aa <- c()
	for (fi in fl) {
		tmp <- read.table(paste0(dir, fi), header = T, stringsAsFactors = F)
		aa <- rbind(aa, tmp[1,])
	}
	aa$Cas <- Cas
	aa$Trait <- names
	res <- rbind(res, aa)
	cat(Cas, '\t')
}


## classify categories
anthro <- c('Body mass index', 'Birth weight', 'Waist circumference', 'Waist-hip ratio', 'Height') #5
mental <- c('Alzheimer\'s disease', 'Bipolar disorder', 'Major depression disorder',
			'Attention deficit hyperactivity disorder', 'Autism spectrum disorder', 'Schizophrenia') #6
metabolic <- c('Fasting glucose', 'Fasting insulin', 'Fasting proinsulin', 'HDL cholesterol',
				'LDL cholesterol', 'Total cholesterol', 'Triglycerides') #7
disease <- c('Coronary artery disease', 'Celiac disease', 'Inflammatory bowel disease',
				'Ulcerative colitis', 'Rheumatoid arthritis', 'Type 2 diabetes', 'Breast cancer') #7
others <- c('Ever smoke', 'Educational attainment', 'Lifespan', 'Telomere length') #4


cate <- c('anthro', 'mental', 'metabolic', 'disease', 'others')

res$category <- NA
for (i in 1:length(cate)) {
	res[res$Trait %in% get(cate[i]),'category'] <- cate[i]
}

saveRDS(res, file = '/opt/working/projects/prj_026_CRISPR/Rev2203/Cas_all_heritability_enrichment_results.rds')
write.table(res, file = '/opt/working/projects/prj_026_CRISPR/Rev2203/Cas_all_heritability_enrichment_results.txt',
				quote = F, sep = '\t', row.names = F)


## output results to TableS3

anthro <- c('Body mass index', 'Birth weight', 'Waist circumference', 'Waist-hip ratio', 'Height') #5
disease <- c('Coronary artery disease', 'Celiac disease', 'Inflammatory bowel disease',
				'Ulcerative colitis', 'Rheumatoid arthritis', 'Type 2 diabetes', 'Breast cancer') #7
mental <- c('Alzheimer\'s disease', 'Bipolar disorder', 'Major depression disorder',
			'Attention deficit hyperactivity disorder', 'Autism spectrum disorder', 'Schizophrenia') #6
metabolic <- c('Fasting glucose', 'Fasting insulin', 'Fasting proinsulin', 'HDL cholesterol',
				'LDL cholesterol', 'Total cholesterol', 'Triglycerides') #7
others <- c('Ever smoke', 'Educational attainment', 'Lifespan', 'Telomere length') #4

trait_new_order <- c(anthro[order(anthro)], disease[order(disease)], mental[order(mental)], metabolic[order(metabolic)], others[order(others)])



dd <- readRDS('~/Documents/Projects/NGG/Revision2203-debug/Cas_all_heritability_enrichment_results.rds')

cas_info <- read.csv('~/Documents/Projects/NGG/Revision-1/CRISPR_PAM_full_table2201.csv', header = T, stringsAsFactors = F)
cas <- unique(cas_info$Cas)

dd$category[dd$category == 'anthro'] <- 'Anthropometric'
dd$category[dd$category == 'disease'] <- 'Disease'
dd$category[dd$category == 'mental'] <- 'Mental disorder'
dd$category[dd$category == 'metabolic'] <- 'Metabolite'
dd$category[dd$category == 'others'] <- 'Others'


dd$Cas <- factor(dd$Cas, levels = cas)
dd$Trait <- factor(dd$Trait, levels = trait_new_order)

ddd <- dd[order(dd$Cas, dd$Trait),]

ddd <- ddd[, c(8:10,2:7)]
colnames(ddd)[3] <- 'Trait Category'

ddd$FDR <- p.adjust(ddd$Enrichment_p, method = 'bonferroni')

write.table(ddd, file = '~/Documents/Projects/NGG/Revision2203-debug/Cas_TableS3_enrichment_res.txt', row.names = F, quote = F, sep = '\t')













