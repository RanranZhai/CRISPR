## this code is running on BigBlack in R ***********************

Args <- commandArgs(trailingOnly = T)

#m = as.numeric(Args[3])
#m = nchar(Args[1]) ## pam length


GG <- read.table(paste0('/mnt2/CRISPR/', Args[1], '_hg19_pos.tsv'), sep = '\t', stringsAsFactors=F)
CC <- read.table(paste0('/mnt2/CRISPR/', Args[2], '_hg19_pos.tsv'), sep = '\t', stringsAsFactors=F)

#GG <- read.table(paste0('/mnt2/CRISPR/','GG', '_hg19_pos.tsv'), sep = '\t', stringsAsFactors=F)
#CC <- read.table(paste0('/mnt2/CRISPR/', 'CC', '_hg19_pos.tsv'), sep = '\t', stringsAsFactors=F)


all <- rbind(GG, CC)

load('/opt/ShaoYa/Users/ranran/PAM/GRC37_genome/CHRs_length.RData')
N <- 20000
M <- dd[1,23]/N
M <- round(M)


xx_GG <- c()
n <- 0
for (i in 1:22) {
	bb <- all[all$V1 == i,]

	if (i > 1) {
		idx <- which(bb$V2 <= n)
		#bb <- bb[-idx,]
		xx_GG[nrow(xx_GG),4] <- xx_GG[nrow(xx_GG),4] + length(idx)
	} 

	x <- data.frame()

	m <- ceiling((dd[1,i]-n)/M) ## Number of segments on this chromosome
	for (j in 1:m) {
		cut_start <- M*(j-1) + 1 + n
		cut_end <- min(((M*j) + n), dd[1,i])

		x[j,1] <- i
		x[j,2] <- cut_start
		x[j,3] <- cut_end
		x[j,4] <- length(which(bb$V2 >= cut_start & bb$V2 <= cut_end))
	}
	xx_GG <- rbind(xx_GG, x)

	n <- (dd[1,i]-n)%%M
	n <- M-n ## left positions from current chromosome

	cat(i, ' chromosome need ', n, ' from the next chromosome', '\t')
}


colnames(xx_GG) <- c('Chr', 'cut_start', 'cut_end', Args[1])
write.table(xx_GG, file = paste0('/opt/working/projects/prj_026_CRISPR/Rev2203/PAM_auto_2203/', Args[1], '_PAM_count_both_strand.txt'), quote = F, sep = '\t', row.names = F)



## output X segments count

#bb <- all[all[,1] == 'X',]
#x <- data.frame()
#m <- ceiling(155270560/M) ## Number of segments on this chromosome, 155270560 is the length of the X chromosome
#for (j in 1:m) {
#	x[j,1] <- nrow(bb[which(bb[,2] > M*(j-1) & bb[,2] < M*j),])
#	x[j,2] <- M*(j-1) + 1 
#	x[j,3] <- min((M*j), dd[1,i])
#}
#x$CHR <- 'chrX'
#xx <- x
#
#
#xx <- xx[,c(4,2,3,1)]
#colnames(xx) <- c('Chr', 'cut_start', 'cut_end', Args[1])
#write.table(xx, file = paste0('/opt/working/projects/prj_026_CRISPR/Rev2203/PAM_X_2203/', Args[1], '_PAM_count_both_strand.txt'), quote = F, sep = '\t', row.names = F)






