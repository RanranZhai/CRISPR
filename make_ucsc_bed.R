

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~make ucsc bed files for ldsc make.annot.py~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## to run this script, 4 arguments will be needed
## Args[1] : file name of input, 
##			which needs three columns 
##			column 1 indicating chromsome (e.g. 1)
##			column 2 start position
##			column 3 end position (optional)
##
##Args[2] : number, indicating bed file start from original start
##
##Args[3] : number, indicating bed file end from original start
##
##Args[4] : file name of output (result)



Args <- commandArgs(trailingOnly = T)

data <- read.table(gzfile(Args[1]), sep = '\t', header = F)


start <- as.numeric(Args[2])
end <- as.numeric(Args[3])

data[,1] <- paste0('chr', data[,1])

data[,2] <- data[,2] + start
data[,3] <- data[,2] + end

write.table(data, file = Args[4], sep = '\t', quote = F, row.names = F, col.names = F)









