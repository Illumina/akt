mendel <- read.table(commandArgs(TRUE)[1],header=T,colClasses=c(rep("character",6),rep("numeric",6)))
## mendel$type1 <-apply(mendel[,c('DAD_GT','MUM_GT')],1, function(x) paste(sort(as.character(x)),collapse="_"))
## mendel$type2 <-apply(mendel[,c('DAD_GT','MUM_GT')],1, function(x) paste(as.character(x),collapse="_"))

mendel$type1 <-apply(mendel[,c('DAD_GT','MUM_GT')],1, function(x) paste(sort(as.character(x)),collapse="_"))
mendel$type2 <-apply(mendel[,c('DAD_GT','MUM_GT')],1, function(x) paste(as.character(x),collapse="_"))


mendel.trio <- subset(mendel,DAD_GT!="."&MUM_GT!=".")
if(NROW(mendel.trio)>0) {
    rr <- tapply(mendel.trio$CHILD_RR,mendel.trio$type1,sum)
    tab <- data.frame(PARENT1=gsub("_.*","",names(rr)),PARENT2=gsub(".*_","",names(rr)),stringsAsFactors=FALSE)
    tab$RR <- rr
    tab$RA <- tapply(mendel.trio$CHILD_RA,mendel.trio$type1,sum)
    tab$AA <- tapply(mendel.trio$CHILD_AA,mendel.trio$type1,sum)
    tab$error_rate <- 100*tapply(mendel.trio$NERROR,mendel.trio$type1,sum)/(tab$RR+tab$RA+tab$AA)
    tab$het_rate <- 100*tab$RA/(tab$RR+tab$RA+tab$AA)
    cat("\nTrio summary:\n")
    print(tab,row.names=FALSE,digits=3)
    cat("\n")
    cat("\n% error excluding RR-RR-RR:",100*sum(tapply(mendel.trio$NERROR,mendel.trio$type1,sum) )/(sum(tab$RR[!(tab$DAD=="RR"&tab$MUM=="RR")])+sum(tab$RA+tab$AA)))
    cat("\n")
    cat("\n")
    
}
    
    
mendel.duo <- subset(mendel,DAD_GT=="."|MUM_GT==".")
if(NROW(mendel.duo)>0) {
    rr <- tapply(mendel.duo$CHILD_RR,mendel.duo$type1,sum)
    tab <- data.frame(DAD=gsub("_.*","",names(rr)),MUM=gsub(".*_","",names(rr)))
    tab$RR <- rr
    tab$RA <- tapply(mendel.duo$CHILD_RA,mendel.duo$type1,sum)
    tab$AA <- tapply(mendel.duo$CHILD_AA,mendel.duo$type1,sum)
    tab$error_rate <- 100*tapply(mendel.duo$NERROR,mendel.duo$type1,sum)/(tab$RR+tab$RA+tab$AA)
    tab$het_rate <- 100*tab$RA/(tab$RR+tab$RA+tab$AA)
    cat("\nDuo summary:\n")
    print(tab,row.names=FALSE)
}

## par(mar=c(4,4,1,1))
## boxplot(mendel.duo$HET_RATE~mendel.duo$type,col='light grey',ylim=c(0,1))
## grid()
## abline(0.5,0,col=2)


