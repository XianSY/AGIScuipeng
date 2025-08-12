#!/opt/conda/envs/emmax-anno/bin/Rscript --verbose
args <- commandArgs(TRUE)

eigenvec_path = args[1]

names<-strsplit(eigenvec_path,"/")
eigenvec_names = names[[1]][length(names[[1]])]
eigenvec_names = strsplit(eigenvec_names,"\\.")[[1]][1]


eigenvec_data = read.table(eigenvec_path,header=F)
eigenvec_data = eigenvec_data[,c(1,2,3,4,5)]

names(eigenvec_data) = c("ID1","ID2","pca1","pca2","pca3")


png( 
  filename = paste(eigenvec_names,'png',sep='.'),
  width = 1000,
  height = 800,
  units = 'px',
  bg = 'white'
)
plot(eigenvec_data$pca1,eigenvec_data$pca2,pch = c(8),col = c("purple"),main="PCA",xlab="pca1",ylab = "pca2")
dev.off()

print(paste(eigenvec_names,'png',sep='.'))


#ggsave(paste(eigenvec_names,'png',sep='.'),p,width=10,height = 10,dpi=300)
#ggsave(paste(eigenvec_names,'pdf',sep='.'),p,width=10,height = 10,dpi=300)
