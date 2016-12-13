counts <- read.table(commandArgs(TRUE)[1],header=T,as.is=TRUE)
n <- sum(rowSums(counts[,c("RR","RA","AA")]))
x <- sum(counts[,"RA"])

print(binom.test(x,n))


