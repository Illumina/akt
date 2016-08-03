cn <- c("id",paste0("PC",1:20),'eval')
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
fname <- commandArgs(TRUE)[1]
pc <- read.table(fname,as.is=T,col.names=head(cn,-1))

pdf(paste0(fname,".pdf"),width=14)
par(mfrow=c(1,2),mar=c(4,4,1,1),cex=1.2)
plot(pc$PC1~pc$PC2,xlab='PC2',ylab='PC1')
plot(pc$PC3~pc$PC2,xlab='PC2',ylab='PC3')
dev.off()

