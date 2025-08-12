#!/opt/conda/envs/emmax-anno/bin/Rscript --verbose
args <- commandArgs(TRUE)

eigenvec_path = args[1]

names<-strsplit(eigenvec_path,"/")
cverror_names = names[[1]][length(names[[1]])]
cverror_names = strsplit(cverror_names,"\\.")[[1]][1]


cverror_data = read.table(eigenvec_path,header=F)


names(cverror_data) = c("k","CV")


png( 
  filename = paste(cverror_names,'png',sep='.'),
  width = 1000,
  height = 800,
  units = 'px',
  bg = 'white'
)
plot(cverror_data$k,cverror_data$CV,,type = "b",pch = c(21),col = c("black"),bg = "black",lwd = 2,font = 2,font.lab = 2,font.axis = 2,xlab="K value",ylab="Cross-Validation (CV) errors")
dev.off()

print(paste(cverror_names,'png',sep='.'))


