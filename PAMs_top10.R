

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









