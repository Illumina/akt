cols <-  c("#E41A1CCC","#377EB8CC","#4DAF4ACC","#984EA3CC","#FF7F00CC")
cn <- c("id",paste("PC",1:20,sep=""),'eval')

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
ogp <- read.table(paste(script.basename,"/../data/","integrated_call_samples_v3.20130502.ALL.panel",sep=""),header=T,colClasses='character')
print(head(ogp))
ogp.pc <- read.table(paste(script.basename,"/../data/wgs.1000G.phase3.pca",sep=""),col.names=head(cn,21))
names(cols) <- unique(ogp$super_pop)
    
fname <- commandArgs(TRUE)[1]
pc <- read.table(fname,as.is=T,col.names=head(cn,-1))

add_text <- nrow(pc)<10

a <- .5
outfile=paste(fname,".pdf",sep="")
pdf(outfile,width=14)
par(mfrow=c(1,2),mar=c(4,4,1,1),cex=1.2)
plot(ogp.pc$PC1~ogp.pc$PC2,pch=16,col=cols[ogp$super_pop],xlab='PC2',ylab='PC1')
points(pc$PC1~pc$PC2,pch=16)
if(add_text)
{
    text(pc$PC2-a,pc$PC1-a,labels=pc[,1],adj=c(1,1))
}
legend("bottomleft",names(cols),pch=16,col=cols)
plot(ogp.pc$PC3~ogp.pc$PC2,pch=16,col=cols[ogp$super_pop],xlab='PC2',ylab='PC3')
points(pc$PC3~pc$PC2,pch=16)
if(add_text)
{
    text(pc$PC2-a,pc$PC3-a,labels=pc[,1],adj=c(1,1))
}
dev.off()

cmd=paste("convert -density 150 ",outfile,paste(fname,".png",sep=""),collapse=" ")
print(cmd)
system(cmd)
       
