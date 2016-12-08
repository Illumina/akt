counts <- read.table(commandArgs(TRUE)[1],header=T,as.is=TRUE)
counts <-     subset(counts,DAD_GT=="RA" | MUM_GT=="RR")
for(i in 0:2)
{
    x <- sum(counts[,paste("RA_",i,sep="")])
    n <- sum(counts[,paste(c("RR_","RA_","AA_"),i,sep="")])
    if(n>0)
    {
        print(paste("status ",i,sep=""))
        print(binom.test(x,n))
    }
}

