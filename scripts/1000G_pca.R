cols <-  c("#E41A1CCC","#377EB8CC","#4DAF4ACC","#984EA3CC","#FF7F00CC")
cn <- c("id",paste0("PC",1:20),'eval')

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
ogp <- read.table(paste0(script.basename,"/../data/","integrated_call_samples_v3.20130502.ALL.panel"),header=T,colClasses='character')
print(head(ogp))
ogp.pc <- read.table(paste0(script.basename,"/../data/1000G.phase3.pca"),col.names=head(cn,21))
names(cols) <- unique(ogp$super_pop)
    
fname <- commandArgs(TRUE)[1]
pc <- read.table(fname,as.is=T,col.names=head(cn,-1))

pdf(paste0(fname,".pdf"),width=14)
par(mfrow=c(1,2),mar=c(4,4,1,1),cex=1.2)
plot(ogp.pc$PC1~ogp.pc$PC2,pch=16,col=cols[ogp$super_pop],xlab='PC2',ylab='PC1')
points(pc$PC1~pc$PC2)
legend("topleft",names(cols),pch=16,col=cols)
plot(ogp.pc$PC3~ogp.pc$PC2,pch=16,col=cols[ogp$super_pop],xlab='PC2',ylab='PC3')
points(pc$PC3~pc$PC2)
dev.off()

