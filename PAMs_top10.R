

Args <- commandArgs(trailingOnly = T)

## get chromosome length
load('/Volumes/ShaoYa/Users/ranran/PAM/GRC37_genome/CHRs_length.RData')

#load('/opt/ShaoYa/Users/ranran/PAM/GRC37_genome/CHRs_length.RData')




N <- 20000
M <- dd[1,23]/N
M <- round(M)


dir.create(paste0('/Volumes/ShaoYa/Users/ranran/PAM/', Args[1], '/top10'))
dir.create(paste0('/Volumes/ShaoYa/Users/ranran/PAM/', Args[1], '/top10/annotation'))
dir.create(paste0('/Volumes/ShaoYa/Users/ranran/PAM/', Args[1], '/top10/results'))



xx <- c()
n <- M
for (i in 1:22) {
	aa <-  read.table(paste0('/Volumes/ShaoYa/Users/ranran/PAM/', Args[1], '/bed/', Args[1], '.', i, '.bed'), header = F, stringsAsFactors = F)
	bb <- aa[order(aa[,2]),]
	rownames(bb) <- c(1:nrow(bb))

	if (i > 1) {
		xx[nrow(xx),1] <- xx[nrow(xx),1] + nrow(bb[which(bb[,2] < (M-n)),])
	} 

	if (length(which(bb[,2] < (M-n))) > 0) {
		bb <- bb[-which(bb[,2] < (M-n)),]
	}
	x <- data.frame()

	m <- ceiling((dd[1,i]-(M-n))/M) ## Number of segments on this chromosome
	for (j in 1:m) {
		x[j,1] <- nrow(bb[which(bb[,2] > M*(j-1) & bb[,2] < M*j),])
		x[j,2] <- M*(j-1) + 1 + (M-n)
		x[j,3] <- min((M*j + (M-n)), dd[1,i])
	}
	x$CHR <- i
	xx <- rbind(xx, x)

	n <- (dd[1,i]-(M-n))%%M ## left positions from current chromosome

	cat(i, '\n')

}

yy <- xx[order(xx[,1], decreasing = T),]
top10 <- yy[1:2000,]


idx <- which(top10[,3]-top10[,2] < M-1)

if (length(idx) > 0) {
	for (i in 1:length(idx)) {
		ccc <- top10[idx[i],3]-top10[idx[i],2]
		tmp <- c(0,1, (M-ccc), (top10[idx[i],4]+1))
		top10 <- rbind(top10, tmp)
	}
}



options(scipen = 200)
Top10 <- top10[,c(4,2,3)]
Top10$CHR <- paste0('chr', Top10$CHR)

write.table(Top10, file = paste0('/Volumes/ShaoYa/Users/ranran/PAM/', Args[1], '/top10/Top10.bed'), col.names = F, row.names = F, quote = F, sep = '\t')



#pdf(paste0('/Volumes/ShaoYa/Users/ranran/PAM/hist/Number_of_', 'TTTA', '_hist.pdf'))
#hist(xx[,1], xlab = paste0('Number of ', 'TTTA'))
#abline(v=min(top10[1:2000,1]), col = 2)
#dev.off()
#
#write.table(Top10, file = paste0('/Volumes/ShaoYa/Users/ranran/PAM/', 'TTTA', '/top10/Top10.bed'), col.names = F, row.names = F, quote = F, sep = '\t')




dir <- list.dirs('/Volumes/ShaoYa/Users/ranran/PAM/')
dir <- dir[grep('results', dir)]


#names <- c('Pleio', 'ALZ', 'BMI', 'BREAST', 'BW3', 'CAD2015', 'EDU2016', 'EUR_CD', 'EUR.IBD', 'EUR.UC', 
#	'FASTGLU', 'FASTINS', 'FASTPROINS', 'HDL', 'Height', 'LDL', 'Lifespan', 'MDD2018', 'PGC_ADHD', 'PGC_ASD', 
#	'RA_GWASmeta', 'SCZ2', 'SMKCPD', 'SMKEVER', 'TC', 'TG', 'TL', 'WC', 'WHR', 'BIP', 'T2D', 'pleio')

names <- c('Pleio', 'Alzheimer\'s disease', 'Body mass index', 'Breast cancer', 'Birth weight', 'Coronary artery disease', 'Educational attainment', 'Celiac disease', 'Inflammatory bowel disease', 'Ulcerative colitis', 
	'Fasting glucose', 'Fasting insulin', 'Fasting proinsulin', 'HDL cholesterol', 'Height', 'LDL cholesterol', 'Lifespan', 'Major depression disorder', 'Attention deficit hyperactivity disorder', 'Autism spectrum disorder', 
	'Rheumatoid arthritis', 'Schizophrenia', 'SMKCPD', 'Ever smoke', 'Total cholesterol', 'Triglycerides', 'Telomere length', 'Waist circumference', 'Waist-hip ratio', 'Bipolar disorder', 'Type 2 diabetes', 'pleio')



ddd <- c()
for (i in 1:(length(dir))) {
	fl <- list.files(dir[i])
	fl <- fl[grep('sumstats.gz.results', fl)]
	dd <- c()
	for (j in 1:length(fl[1:32])) {
		tmp <- read.table(paste0(dir[i], '/', fl[j]), header = T, stringsAsFactors = F)
		dd <- rbind(dd, tmp[1,])
	}
	dd$Trait <- names
	pam <- unlist(strsplit(fl[1], 'V2_'))[2]
	pam <- unlist(strsplit(pam, '_2_formatted.'))[1]
	#pam <- gsub('/', '_', pam)
	dd$PAM <- rep(pam, 32)
	ddd <- rbind(ddd, dd)
	cat(pam, '\n')
}


nejm.colors <- function (n = 8, alpha = 1){
    pal <- c("#BC3C29", "#0072B5", "#E18727", "#20854E", "#7876B1",
             "#6F99AD", "#FFDC91", "#EE4C97")
    acode <- substr(rgb(0, 0, 0, alpha = alpha), 8, 9)
    return(paste0(pal, acode)[1:n])
}

idx <- which(ddd$PAM %in% c('GAW', 'TTTV'))

#colors <- c(nejm.colors(8), nejm.colors(1, .8))

require(ggplot2)

#pdf('/Volumes/ShaoYa/Users/ranran/PAM/heritability_enrichment_all.pdf', 12, 8)
#ggplot(ddd, aes(x = Trait, y = Enrichment, fill = PAM)) + geom_bar(stat = 'identity', position = 'dodge') +
#	geom_linerange(aes(ymin = Enrichment - 1.96*Enrichment_std_error, ymax = Enrichment + 1.96*Enrichment_std_error), position = position_dodge(.9)) +
#	theme(axis.text.x = element_text(angle = 45, vjust = .5)) +  
#	scale_fill_manual(values = colors) +
#	geom_hline(yintercept = 1, col = nejm.colors(1), linetype = 'dashed')
#dev.off()



### split traits into different category



names <- c('Pleio', 'ALZ', 'BMI', 'BREAST', 'BW3', 'CAD2015', 'EDU2016', 'EUR_CD', 'EUR.IBD', 'EUR.UC', 
	'FASTGLU', 'FASTINS', 'FASTPROINS', 'HDL', 'Height', 'LDL', 'Lifespan', 'MDD2018', 'PGC_ADHD', 'PGC_ASD', 
	'RA_GWASmeta', 'SCZ2', 'SMKCPD', 'SMKEVER', 'TC', 'TG', 'TL', 'WC', 'WHR', 'BIP', 'T2D', 'pleio')

names <- c('Pleio', 'Alzheimer\'s disease', 'Body mass index', 'Breast cancer', 'Birth weight', 'Coronary artery disease', 'Educational attainment', 'Celiac disease', 'Inflammatory bowel disease', 'Ulcerative colitis', 
	'Fasting glucose', 'Fasting insulin', 'Fasting proinsulin', 'HDL cholesterol', 'Height', 'LDL cholesterol', 'Lifespan', 'Major depression disorder', 'Attention deficit hyperactivity disorder', 'Autism spectrum disorder', 
	'Rheumatoid arthritis', 'Schizophrenia', 'SMKCPD', 'Ever smoke', 'Total cholesterol', 'Triglycerides', 'Telomere length', 'Waist circumference', 'Waist-hip ratio', 'Bipolar disorder', 'Type 2 diabetes', 'pleio')


anthro <- c('Body mass index', 'Birth weight', 'Waist circumference', 'Waist-hip ratio', 'Height') #5

mental <- c('Alzheimer\'s disease', 'Bipolar disorder', 'Major depression disorder', 'Attention deficit hyperactivity disorder', 'Autism spectrum disorder', 'Schizophrenia') #6

metabolic <- c('Fasting glucose', 'Fasting insulin', 'Fasting proinsulin', 'HDL cholesterol', 'LDL cholesterol', 'Total cholesterol', 'Triglycerides') #7

disease <- c('Coronary artery disease', 'Celiac disease', 'Inflammatory bowel disease', 'Ulcerative colitis', 'Rheumatoid arthritis', 'Type 2 diabetes', 'Breast cancer') #7

others <- c('Pleio', 'pleio', 'SMKCPD', 'Ever smoke', 'Educational attainment', 'Lifespan', 'Telomere length') #7


cate <- c('anthro', 'mental', 'metabolic', 'disease', 'others')

ddd$category <- NA
for (i in 1:length(cate)) {
	ddd[ddd$Trait %in% get(cate[i]),'category'] <- cate[i]
}

save(ddd, file = '/Volumes/ShaoYa/Users/ranran/PAM/PAM_heritability_enrichment.RData')



load('~/Documents/Projects/NGG/PAM_heritability_enrichment.RData')
idx <- which(ddd$Trait %in% c('Fasting proinsulin', 'Pleio', 'pleio', 'SMKCPD'))
#idx2 <- which(ddd$PAM %in% c('GAW_top10', 'TTTV_top10'))
idx2 <- which(ddd$PAM %in% c('GAW_top10', 'TTTV_top10', 'GC_content'))

dd <- ddd[-c(idx, idx2),]

nejm.colors <- function (n = 8, alpha = 1){
    pal <- c("#BC3C29", "#0072B5", "#E18727", "#20854E", "#7876B1",
             "#6F99AD", "#FFDC91", "#EE4C97")
    acode <- substr(rgb(0, 0, 0, alpha = alpha), 8, 9)
    return(paste0(pal, acode)[1:n])
}

colors <- c(nejm.colors(8), nejm.colors(8, .8), nejm.colors(8, .6), nejm.colors(5, .4))


idx <- grep('top10', dd$PAM)
dd$PAM[idx] <- unlist(strsplit(dd$PAM[idx], '_top10'))

dd$PAM[which(dd$PAM == 'GAG')] <- 'NGAG'

require(ggplot2)
#pdf('~/Documents/Projects/NGG/PAM_heritability_enrichment_20210111.pdf', 7,8)
pdf('~/Documents/Projects/NGG/PAM_heritability_enrichment_20210209.pdf', 7,8)
ggplot(dd, aes(x = PAM, y = Enrichment, fill = Trait)) + geom_bar(stat = 'identity', position = 'dodge') +
	geom_hline(yintercept = 1, col = nejm.colors(1), linetype = 'dashed') +
	geom_linerange(aes(ymin = Enrichment - 1.96*Enrichment_std_error, ymax = Enrichment + 1.96*Enrichment_std_error), position = position_dodge(.9), size = .3) +
	theme(axis.text.x = element_text(angle = 45, vjust = .5)) + 
	theme_classic() + 
	scale_fill_manual(values = colors) +
	#geom_hline(yintercept = 1, col = nejm.colors(1), linetype = 'dashed') +
	facet_grid(rows = vars(category)) + 
	theme(legend.position = 'none') + 
	theme(strip.text.y = element_blank(), axis.title.x = element_blank())
dev.off()


require(ggplot2)
pdf('/Volumes/ShaoYa/Users/ranran/PAM/PAM_heritability_enrichment_20210119.pdf', 7,8)
ggplot(dd, aes(x = PAM, y = Enrichment, fill = Trait)) + geom_bar(stat = 'identity', position = 'dodge') +
	geom_hline(yintercept = 1, col = nejm.colors(1), linetype = 'dashed') +
	geom_linerange(aes(ymin = Enrichment - 1.96*Enrichment_std_error, ymax = Enrichment + 1.96*Enrichment_std_error), position = position_dodge(.9), size = .3) +
	theme(axis.text.x = element_text(angle = 45, vjust = .5)) + 
	theme_classic() + 
	scale_fill_manual(values = colors) +
	#geom_hline(yintercept = 1, col = nejm.colors(1), linetype = 'dashed') +
	facet_grid(rows = vars(category)) + 
	theme(legend.position = 'none') + 
	theme(strip.text.y = element_blank(), axis.title.x = element_blank())
dev.off()





traits <- dd$Trait[!duplicated(dd$Trait)]
traits <- traits[order(traits)]


cate <- c('anthro', 'mental', 'metabolic', 'disease', 'others')
Cate <- c('Anthropometric', 'Mental disease', 'Metabolite', 'Disease', 'Others')
for (i in 1:length(cate)) {
	ddd <- dd[which(dd$category == cate[i]),]
	idx <- which(traits %in% ddd$Trait)
	pdf(paste0('~/Documents/Projects/NGG/PAM_heritability_p', i, '_legend.pdf'))
	p <- ggplot(ddd, aes(x = PAM, y = Enrichment, fill = Trait)) + geom_bar(stat = 'identity', position = 'dodge') +
	geom_linerange(aes(ymin = Enrichment - 1.96*Enrichment_std_error, ymax = Enrichment + 1.96*Enrichment_std_error), position = position_dodge(.9)) +
	theme(axis.text.x = element_text(angle = 45, vjust = .5)) + 
	theme_classic() + 
	scale_fill_manual(values = colors[idx]) + guides(fill = guide_legend(title = Cate[i]))
	print(p)
	dev.off()
}



### venn plot for different CRISPR-Cas proteins
require(VennDiagram)
require(RColorBrewer)

load('~/Documents/Projects/NGG/annotated_seg.RData')

load('~/Documents/Projects/NGG/GAAG_GAGA/annotated_seg.RData')


#myCol <- brewer.pal(3, "Pastel2")


## Cas9 (NGCG, NGAG, NGG, and GC content)
set1 <- cc[cc$PAM == 'GC_content', 'Start']
set2 <- cc[cc$PAM == 'NGG', 'Start']
set3 <- cc[cc$PAM == 'GAG', 'Start']
set4 <- cc[cc$PAM == 'NGCG', 'Start']

myCol <- brewer.pal(4, "Pastel2")

venn.diagram(
    x = list(set1, set2, set3,set4),
    category.names = c("GC_content", "NGG", "NGAG", "NGCG"),
    filename = '~/Documents/Projects/NGG/venn_plot/SpCas9_and_var.png',
    output=TRUE, fill = myCol
)


## xCas9 (GAA, GAT)
set1 <- cc[cc$PAM == 'GAA', 'Start']
set2 <- cc[cc$PAM == 'GAT', 'Start']

myCol <- brewer.pal(3, "Pastel2")

venn.diagram(
    x = list(set1, set2),
    category.names = c('GAA', 'GAT'),
    filename = '~/Documents/Projects/NGG/venn_plot/xCas9.png',
    output=TRUE, fill = myCol[1:2]
)


## xCas9 (GAA, GAT)
set1 <- cc[cc$PAM == 'TTTA', 'Start']
set2 <- cc[cc$PAM == 'TTTC', 'Start']
set3 <- cc[cc$PAM == 'TTTG', 'Start']

myCol <- brewer.pal(3, "Pastel2")

venn.diagram(
    x = list(set1, set2, set3),
    category.names = c('TTTA', 'TTTC', 'TTTG'),
    filename = '~/Documents/Projects/NGG/venn_plot/cpf1.png',
    output=TRUE, fill = myCol
)



set1 <- cc[cc$PAM == 'GC_content', 'Start']
set2 <- cc[cc$PAM == 'GAGA', 'Start']
set3 <- cc[cc$PAM == 'GAAG', 'Start']

myCol <- brewer.pal(3, "Pastel2")

venn.diagram(
    x = list(set1, set2, set3),
    category.names = c("GC_content", "GAGA", "GAAG"),
    filename = '~/Documents/Projects/NGG/venn_plot/GAAG_GAGA.png',
    output=TRUE, fill = myCol
)


set1 <- cc[cc$PAM == 'GC_content', 'Start']
set2 <- cc[cc$PAM == 'GTGT', 'Start']
set3 <- cc[cc$PAM == 'GTTG', 'Start']

myCol <- brewer.pal(3, "Pastel2")

venn.diagram(
    x = list(set1, set2, set3),
    category.names = c("GC_content", "GTGT", "GTTG"),
    filename = '~/Documents/Projects/NGG/venn_plot/GTTG_GTGT.png',
    output=TRUE, fill = myCol
)



set1 <- cc[cc$PAM == 'GC_content', 'Start']
set2 <- cc[cc$PAM == 'GCCG', 'Start']
set3 <- cc[cc$PAM == 'GCGC', 'Start']

myCol <- brewer.pal(3, "Pastel2")

venn.diagram(
    x = list(set1, set2, set3),
    category.names = c("GC_content", "GCGC", "GCCG"),
    filename = '~/Documents/Projects/NGG/venn_plot/GCCG_GCGC.png',
    output=TRUE, fill = myCol
)



set1 <- cc[cc$PAM == 'GC_content', 'Start']
set2 <- cc[cc$PAM == 'GAGA', 'Start']
set3 <- cc[cc$PAM == 'GCGC', 'Start']
set4 <- cc[cc$PAM == 'GTGT', 'Start']


myCol <- brewer.pal(4, "Pastel2")

venn.diagram(
    x = list(set1, set2, set3, set4),
    category.names = c("GC_content", 'GAGA', "GCGC", "GTGT"),
    filename = '~/Documents/Projects/NGG/venn_plot/GNGN.png',
    output=TRUE, fill = myCol
)



set1 <- cc[cc$PAM == 'GC_content', 'Start']
set2 <- cc[cc$PAM == 'GGGG', 'Start']
set3 <- cc[cc$PAM == 'GCGC', 'Start']
set4 <- cc[cc$PAM == 'GCCG', 'Start']


myCol <- brewer.pal(4, "Pastel2")

venn.diagram(
    x = list(set1, set2, set3, set4),
    category.names = c("GC_content", 'GGGG', "GCGC", "GCCG"),
    filename = '~/Documents/Projects/NGG/venn_plot/GC.png',
    output=TRUE, fill = myCol
)






### try cowplot


## load('~/Documents/Projects/NGG/PAM_heritability_enrichment.RData')
## idx <- which(ddd$Trait %in% c('FASTPROINS', 'Pleio', 'pleio'))
## idx2 <- which(ddd$PAM %in% c('GAW_top10', 'TTTV_top10'))
## dd <- ddd[-c(idx, idx2),]
## 
## nejm.colors <- function (n = 8, alpha = 1){
##     pal <- c("#BC3C29", "#0072B5", "#E18727", "#20854E", "#7876B1",
##              "#6F99AD", "#FFDC91", "#EE4C97")
##     acode <- substr(rgb(0, 0, 0, alpha = alpha), 8, 9)
##     return(paste0(pal, acode)[1:n])
## }
## 
## 
## require(cowplot)
## 
## 
## cate <- c('anthro', 'mental', 'metabolic', 'disease', 'others')
## n <- 1
## for (i in 1:length(cate)) {
## 	ddd <- dd[which(dd$category == cate[i]),]
## 	head(ddd)
## 	m <- length(ddd$Trait[!duplicated(ddd$Trait)])
## 	cat(m, '\n')
## 	assign(paste0('p', i), ggplot(ddd, aes(x = PAM, y = Enrichment, fill = Trait)) + geom_bar(stat = 'identity', position = 'dodge') +
## 	geom_linerange(aes(ymin = Enrichment - 1.96*Enrichment_std_error, ymax = Enrichment + 1.96*Enrichment_std_error), position = position_dodge(.9)) +
## 	theme(axis.text.x = element_text(angle = 45, vjust = .5)) + 
## 	theme_classic() + 
## 	scale_fill_manual(values = colors[n:(m+n-1)]) +
## 	geom_hline(yintercept = 1, col = nejm.colors(1), linetype = 'dashed'))
## 	cat(n:(m+n-1), '\n')
## 	n <- m+n
## }
## 
## 
## n <- 1
## for (i in 1:length(cate)) {
## 	ddd <- dd[which(dd$category == cate[i]),]
## 	head(ddd)
## 	m <- length(ddd$Trait[!duplicated(ddd$Trait)])
## 	cat(m, '\n')
## 	if (i < 5) {
## 		assign(paste0('p', i), ggplot(ddd, aes(x = PAM, y = Enrichment, fill = Trait)) + geom_bar(stat = 'identity', position = 'dodge') +
## 		geom_linerange(aes(ymin = Enrichment - 1.96*Enrichment_std_error, ymax = Enrichment + 1.96*Enrichment_std_error), position = position_dodge(.9)) +
## 		theme(axis.text.x = element_text(angle = 45, vjust = .5)) + 
## 		theme_classic() + 
## 		scale_fill_manual(values = colors[n:(m+n-1)]) +
## 		geom_hline(yintercept = 1, col = nejm.colors(1), linetype = 'dashed') + 
## 		theme(legend.position = 'none', axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks.x=element_blank()) + ylim(c(-1,5)))
## 	} else {
## 		assign(paste0('p', i), ggplot(ddd, aes(x = PAM, y = Enrichment, fill = Trait)) + geom_bar(stat = 'identity', position = 'dodge') +
## 		geom_linerange(aes(ymin = Enrichment - 1.96*Enrichment_std_error, ymax = Enrichment + 1.96*Enrichment_std_error), position = position_dodge(.9)) +
## 		theme(axis.text.x = element_text(angle = 45, vjust = .5)) + 
## 		theme_classic() + 
## 		scale_fill_manual(values = colors[n:(m+n-1)]) +
## 		geom_hline(yintercept = 1, col = nejm.colors(1), linetype = 'dashed') + 
## 		theme(legend.position = 'none', axis.title.y = element_blank()) + ylim(c(-1,5)))
## 	}
## 	cat(n:(m+n-1), '\n')
## 	n <- m+n
## }
## 
## pdf('~/Documents/Projects/NGG/PAM_heritability_enrichment_cowplot.pdf', 7,8)
## plot_grid(p1,p2,p3,p4,p5, ncol = 1)
## dev.off()
























