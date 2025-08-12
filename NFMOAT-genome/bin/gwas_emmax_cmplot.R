#!/opt/conda/envs/emmax-anno/bin/Rscript --verbose
library("CMplot")
library("stringr")

args <- commandArgs(TRUE)
assocresultfile <- args[1]
names<-strsplit(assocresultfile,"/")
name<-names[[1]][length(names[[1]])]

assocresult<-read.table(assocresultfile,header=F)
assocresult[,c(5,6)]<-as.data.frame(str_split_fixed(assocresult[,1],":",2))

assocplot<-assocresult[c(1,5,6,4)]
colnames(assocplot)<-c("ID","chr","position","P")
thresholds <- c(0.01,0.05)/nrow(assocplot)

sigresult <- assocplot[which(assocplot[,"P"]<thresholds[2]),]
write.table(sigresult,paste0(assocresultfile,"_significant.lst"),sep="\t",row.names=F,quote = FALSE,col.names=T)

CMplot(assocplot,
       plot.type = "m",
       #chr.labels=paste("Chr",c(1:22,"X","Y"),sep=""),
       file="png",
       file.name = paste0(name,".manh"),
       dpi=300,
       file.output=TRUE,verbose=TRUE,width=14,height=6,
       threshold =c(thresholds,0.000001), #阈值设定，Bonfferoni法校正=0.01/nrow(Pmap)，通常以1/SNP总数为阈值。
       threshold.col=c('grey','black','black'),#控制阈值线的颜色，注意要对应
       threshold.lty = c(1,2,2),#阈值线的类型（7为实线，5为虚线）
       threshold.lwd = c(1,1,1),#阈值线的宽度，也可以调整对角线的宽度，2最合适）
       amplify = T,#放大显示~显著性的SNP点
       signal.cex = c(1,1,1),#设置有效SNP位点的大小,建议选用2/3号为宜
       signal.pch = c(20,20,20),#设置有效SNP位点的形状，16和19为圆形。
       signal.col =c("red","darkgreen","skyblue"),#有效SNP位点的颜色
       chr.den.col=c("darkgreen","yellow","red")#【添加SNP密度】！！！,若不想添加只需要到前一步就可结束
)

CMplot(assocplot,
       plot.type="q",
       threshold = 0.05,
       conf.int.col=NULL,#绘制QQ图中置信区间的颜色，可以用字符或向量
       box=TRUE,#是否添加边框
       file="png",
       file.name = paste0(name,".QQ"),
       dpi=300,
       file.output=TRUE,
       verbose=TRUE,width=10,height=10
)
