#!/usr/bin/env Rscript
library(optparse)
option_list = list(
  make_option(
    c("--expression","-E"),
    default = NULL,
    metavar = "path",
    type = "character",
    help = "genes expression information"
  ),
  make_option(
    c("--gene_info"),
    default = NULL,
    metavar = "path",
    type = "character",
    help = "gene information of beds"
  )
)

opt_parse = OptionParser(option_list = option_list)
opt = parse_args(opt_parse)

expression_data = read.csv(opt$expression,header = T,row.names = 1)


expression_data_filter = expression_data[rowMeans(expression_data) >0.1,]
rownames_df = colnames(expression_data_filter)
#expression_data_filter = expression_data_filter[rowSums(expression_data_filter==0) < 0.95,]
expression_data_filter = t(apply(expression_data_filter,1,function(x) qqnorm(x,plot.it = F)$x))

colnames(expression_data_filter) = rownames_df

gene_info = read.table(opt$gene_info,header = F)
names(gene_info) = c("gene_ID","CHR","start","end")

expression_data_filter = expression_data_filter[rownames(expression_data_filter) %in% gene_info$gene_ID,]

write.table(expression_data_filter,"expression_filter.tab",quote = F,row.names = T,col.names = T,sep = "\t")

gene_info = gene_info[gene_info$gene_ID %in% rownames(expression_data_filter),]
write.table(gene_info,"gene_info.txt",quote = F,row.names = F,col.names = T,sep = "\t")
