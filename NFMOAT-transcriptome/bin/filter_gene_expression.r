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

counts = read.table(opt$gene_counts,sep = ",",row.names = 1,header = T)

library_size <- colSums(counts)
cpm <- sweep(
    counts,
    2,
    library_size,
    FUN="/"
) * 1000000

keep <- rowSums(cpm > 1) >= ncol(cpm)*0.1


filter_counts = counts[keep,]

zero_filter = filter_counts

zero_filter$Geneid = row.names(zero_filter)
zero_filter = zero_filter[,c(length(zero_filter),1:(length(zero_filter)-1))]

write.table(zero_filter,"filter_gene_expression.txt",quote = F,sep = ",",col.names = T,row.names=F)
