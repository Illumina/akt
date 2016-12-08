counts <- read.table(commandArgs(TRUE)[1],header=T,as.is=TRUE)
counts <-     subset(counts,DAD_GT=="RA" | MUM_GT=="RR")
for(i in 0:2)
{
    x <- sum(counts[,paste0("RA_",i)])
    n <- sum(counts[,paste0(c("RR_","RA_","AA_"),i)])
    if(n>0)
    {
        print(paste("status ",i,sep=""))
        print(binom.test(x,n))
    }
}

