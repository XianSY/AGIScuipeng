#!/usr/bin/env Rscript

#source("http://zzlab.net/GAPIT/GAPIT.library.R")
#source("http://zzlab.net/GAPIT/gapit_functions.txt")

source("/opt/gapit/gapit_functions.txt")
source("/opt/gapit/GAPIT.library.R")

library(optparse)


option_list = list(
  make_option(
    c("--expression","-E"),
    type = "character",
    default = NULL,
    metavar = "path",
    help = "input sample expression"
  ),
  make_option(
    c("--gene_GM","-G"),
    type = "character",
    default = NULL,
    metavar = "path",
    help = "gene information"
  ),
  make_option(
    c("--phenotype","-p"),
    type = "character",
    default = NULL,
    metavar = "path",
    help = "phenotype file"
  ),
  make_option(
    c("--output","-o"),
    type = "character",
    default = "./",
    metavar = "path",
    help = "output path"
  )
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

#读取数据
PFKM = read.table(opt$expression,sep = ",",header = T,row.names = 1)
#row.names(PFKM) = PFKM[,1]
#PFKM = PFKM[,-c(1)]
PFKM = PFKM[rowMeans(PFKM)>0,]
PFKM = PFKM[rowSums(PFKM==0)<0.95,]
PFKM_t = data.frame(t(PFKM))
#colnames(PFKM_t) = row.names(PFKM)
#PFKM_t = PFKM_t[-1,]
PFKM_num = data.frame(lapply(PFKM_t,as.numeric))
PFKM_num$taxa = row.names(PFKM_t)
PFKM_num = PFKM_num[,c(ncol(PFKM_num),1:ncol(PFKM_num)-1)]
PFKM$geneid = rownames(PFKM)
PFKM = PFKM[,c(length(PFKM),1:length(PFKM)-1)]


#基因数据加载
myGM = read.table(opt$gene_GM,header = T,sep=",")
myGM$geneid = rownames(myGM)
myGM=myGM[,c(length(myGM),1:length(myGM)-1)]
names(myGM) = c("geneid","chrom","start")
unique(myGM$chrom)

myGM = merge(myGM,PFKM,by = "geneid")[,c(1,2,3)]


#将TPM值缩放到0-2区间
Quantile<- apply(PFKM_num[,-1],2,  # 2 indicates it is for column and 1 indicates it is for row
                 function(A){min_x=as.numeric(quantile(A,0.05)); #quantile函数是指取分位数。A为counts第二列以后的数据，0.05是指取取第二列第0.05的数
                 max_x=as.numeric(quantile(A,0.95));
                 out<-(2*(A-min_x)/(max_x-min_x));
                 out[out>2]<-2;out[out< 0]<- 0;return(out)})
Quantile.t <- as.data.frame(Quantile)
Quantile.t$gene_symbol = PFKM_num[,1]
myGD <-  Quantile.t[,c(ncol(Quantile.t),1: (ncol(Quantile.t)-1))]



phenotype = read.table(opt$phenotype,sep = ",",header = T)

for (i in c(2:ncol(phenotype))){
	trait_name = names(phenotype)[i]
	mypheno = phenotype[,c(1,i)]
	
	print(dim(mypheno))
	print(dim(myGD))
	myGAPIT <- GAPIT(Y=mypheno, #表型
                 GD=myGD, #基因表达量归一化到0-2
                 GM=myGM, #基因信息 位置和染色体
                 #CV=myCV, #样本协变量信息
                 PCA.total=3,
                 NJtree.group = 3,
                 NJtree.type = 'fan', 
                 model= c("GLM","MLM","FarmCPU","CMLM"),
                 SNP.MAF=0,
                 file.output=T,
		 Multiple_analysis=TRUE
	)

	values = data.frame(myGAPIT$GWAS)
	values$FDR = p.adjust(values$P.value,method='BH')
	write.csv(values,paste0("TWAS_GAPIT.CMLM_",trait_name,".csv"),row.names = F)

}

