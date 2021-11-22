## this code is running on BigBlack in R ***********************

Args <- commandArgs(trailingOnly = T)

#m = as.numeric(Args[3])
m = nchar(Args[1]) ## pam length


GG <- read.table(paste0('/mnt2/CRISPR/', Args[1], '_hg19_pos.tsv'), sep = '\t', stringsAsFactors=F)
CC <- read.table(paste0('/mnt2/CRISPR/', Args[2], '_hg19_pos.tsv'), sep = '\t', stringsAsFactors=F)


new_GG <- GG
new_CC <- CC

new_GG[,1] <- paste0('chr', new_GG[,1])

new_GG[,3] <- new_GG[,2] + m - 1



new_CC[,1] <- paste0('chr', new_CC[,1])

#new_CC[,2] <- new_CC[,3]

new_CC[,3] <- new_CC[,2] + m - 1


load('/opt/ShaoYa/Users/ranran/PAM/GRC37_genome/CHRs_length.RData')

N <- 20000
M <- dd[1,23]/N
M <- round(M)


xx_GG <- c()
n <- M
for (i in 1:22) {
	bb <- GG[GG[,1] == i,]

	rownames(bb) <- c(1:nrow(bb))

	if (i > 1) {
		xx_GG[nrow(xx_GG),1] <- xx_GG[nrow(xx_GG),1] + nrow(xx_GG[which(xx_GG[,2] < (M-n)),])
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
	xx_GG <- rbind(xx_GG, x)

	n <- (dd[1,i]-(M-n))%%M ## left positions from current chromosome

	cat(i, '\t')
}

xx <- xx_GG


## PAM count for the reverse strand
xx_GG <- c()
n <- M
for (i in 1:22) {
	bb <- CC[CC[,1] == i,]

	rownames(bb) <- c(1:nrow(bb))

	if (i > 1) {
		xx_GG[nrow(xx_GG),1] <- xx_GG[nrow(xx_GG),1] + nrow(xx_GG[which(xx_GG[,2] < (M-n)),])
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
	xx_GG <- rbind(xx_GG, x)

	n <- (dd[1,i]-(M-n))%%M ## left positions from current chromosome

	cat(i, '\t')

}


xx[,1] <- xx[,1] + xx_GG[,1]
xx <- xx[,c(4,2,3,1)]
colnames(xx) <- c('Chr', 'cut_start', 'cut_end', Args[1])
write.table(xx, file = paste0('/opt/working/projects/prj_026_CRISPR/PAM/', Args[1], '_PAM_count_both_strand.txt'), quote = F, sep = '\t', row.names = F)










