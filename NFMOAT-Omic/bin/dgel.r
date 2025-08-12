#!/usr/bin/env Rscript
library(optparse)

option_list = list(
  make_option(
    c("--expression","-E"),
    default = NULL,
    type = "character",
    metavar = "path",
    help = "expression files"
  )
)

opt_parser = OptionParser(option_list = option_list)
opt        = parse_args(opt_parser)

if(is.null(opt$expression)){
  print_help(opt_parser)
  stop("Please provide a expression file .",call. = FALSE)
}

expression_data = read.table(opt$expression,header = T,row.names = 1,sep=",")

expression_data_t = data.frame(t(expression_data))

for(i in c(1:ncol(expression_data_t))){
  gene_name = colnames(expression_data_t)[i]
  gene_expression = data.frame(row.names(expression_data_t),row.names(expression_data_t),expression_data_t[,gene_name])
  write.table(gene_expression,file = paste0(gene_name,".tx"),quote = FALSE,row.names = FALSE,col.names = FALSE)
}
