#!/usr/bin/env Rscript


library(FactoMineR)
library(ggplot2)
library(optparse)

option_list <- list(
  make_option(
    c("-i", "--counts_file"),
    type = "character",
    default = NULL    ,
    metavar = "path"   ,
    help = "Counts matrix where rows are Sample and columns are length or sample of counts."
  ),
  make_option(
    c("-g", "--group_file"),
    type = "character",
    default = NULL    ,
    metavar = "path"   ,
    help = "group matrix where rows are Sample and columns are group type."
  ),
  make_option(
    c("-o", "--outdir"),
    type = "character",
    default = './'    ,
    metavar = "path"   ,
    help = "Output directory."
  )
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)
featureCounts = read.delim(opt$counts_file,row.names = 1,sep = "")
gene_info = featureCounts[,c(1:2)]
gene_info$chr_fix = sapply(gene_info$Chr,function(x){strsplit(x,";")[[1]][1]})
gene_info$Start_fix = sapply(gene_info$Start,function(x){strsplit(x,";")[[1]][1]})
write.table(file = "gene_info.txt",sep = ",",gene_info[,c("chr_fix","Start_fix")],row.names = F,quote = FALSE)
featureCounts = featureCounts[,-c(1:4)]
RKB = featureCounts$Length / 1000
featureCounts = featureCounts[,2:length(featureCounts)]/RKB
featureCounts = featureCounts/colSums(featureCounts)*1e6
featureCounts = featureCounts[rowMeans(featureCounts)>0.5,]
featureCounts = featureCounts[rowSums(featureCounts==0)<0.9*ncol(featureCounts),]

featureCounts$geneid = row.names(featureCounts)
featureCounts = featureCounts[,c(length(featureCounts),1:length(featureCounts)-1)]
write.table(file = paste(opt$outdir,"TPM.csv",sep = ""),quote = FALSE,featureCounts,row.name = F,sep = ",")
group = read.csv(opt$group_file,row.names = 1)
gene = t(featureCounts[,-1])
gene.pca = FactoMineR::PCA(gene,ncp = 10,scale.unit = TRUE,graph = FALSE)
plot(gene.pca)
pca = data.frame(gene.pca$ind$coord[,1:2])
pca_eig1 = round(gene.pca$eig[1,2],2)
pca_eig2 = round(gene.pca$eig[2,2],2)
pca$sample = row.names(pca)
group = group[row.names(pca),]
pca = cbind(pca,group)
pca_pic = ggplot2::ggplot(data = pca, aes(x = Dim.1, y = Dim.2, color = group)) + geom_point() +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  xlab(paste("PCA1 ", pca_eig1, "%", sep = "")) + ylab(paste("PCA2 ", pca_eig2, "%")) +
  geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text(aes(label = sample), color = "black", size = 3)
ggsave("PCA.png",pca_pic,width=7,height=7)
