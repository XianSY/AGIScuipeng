#!/usr/bin/env Rscript
library(optparse)

option_list = list(
  make_option(
    c("--gene_counts","-E"),
    default = NULL,
    type = "character",
    metavar = "path",
    help = "expression files"
  )
)


opt_parser = OptionParser(option_list = option_list)
opt        = parse_args(opt_parser)

counts = read.table(opt$gene_counts,sep = " ",row.names = 1,header = T)

filter_counts = counts[rowSums(counts)>1,]

zero_filter = filter_counts[rowSums(filter_counts==0) < ncol(filter_counts)*0.95,]

zero_filter$Geneid = row.names(zero_filter)
zero_filter = zero_filter[,c(length(zero_filter),1:(length(zero_filter)-1))]

write.table(zero_filter,"filter_gene_expression.txt",quote = F,sep = ",",col.names = T,row.names=F)
