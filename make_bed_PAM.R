## this code is running on BigBlack in R ***********************

Args <- commandArgs(trailingOnly = T)



GG <- read.table(paste0('/Volumes/ShaoYa/Users/ranran/Data/PAM_bed/', Args[1], '_hg19_pos.tsv'), sep = '\t', stringsAsFactors=F)
CC <- read.table(paste0('/Volumes/ShaoYa/Users/ranran/Data/PAM_bed/', Args[2], '_hg19_pos.tsv'), sep = '\t', stringsAsFactors=F)


new_GG <- GG
new_CC <- CC

new_GG[,1] <- paste0('chr', new_GG[,1])

new_GG[,3] <- new_GG[,2] - 1

new_GG[,2] <- new_GG[,3] - 20



new_CC[,1] <- paste0('chr', new_CC[,1])

new_CC[,2] <- new_CC[,3] + 1

new_CC[,3] <- new_CC[,2] + 20


dir.create(paste0('/Volumes/ShaoYa/Users/ranran/PAM/', Args[1]))
dir.create(paste0('/Volumes/ShaoYa/Users/ranran/PAM/', Args[1], '/bed'))
options(scipen = 200)
for (i in 1:22) {
	data1 <- new_GG[which(new_GG[,1] == paste0('chr', i)),]
	data2 <- new_CC[which(new_CC[,1] == paste0('chr', i)),]
	write.table(rbind(data1, data2), file = paste0('/Volumes/ShaoYa/Users/ranran/PAM/', Args[1],  '/bed/', Args[1], '.', i, '.bed'), sep = '\t', quote = F, row.names = F, col.names = F)
	cat(i, '\n')
}
options(scipen = 0)








