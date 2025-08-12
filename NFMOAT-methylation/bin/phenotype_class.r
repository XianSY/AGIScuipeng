#!/usr/bin/env Rscript
library(optparse)
option_list = list(
  make_option(
    c("-i","--input_phenotype"),
    type = "character",
    default = NULL,
    metavar = "path",
    help = "phenotype csv file"),
  make_option(
    c("-o","--outdir"),
    type = "character",
    default = './',
    metavar = "path",
    help = "Output phenotype file"
  )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

phenotype = read.csv(opt$input_phenotype,header = TRUE)
for(i in names(phenotype)[-1]){
  new_phenotype = phenotype[,c("sample","sample",i)]
  write.table(new_phenotype,file = paste("phenotype_",i,".txt",sep = ""),quote = FALSE,col.names = TRUE,row.names = FALSE)
}
