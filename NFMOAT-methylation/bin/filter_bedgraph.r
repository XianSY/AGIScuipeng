#!/usr/bin/env Rscript
library(optparse)
library(data.table)

option_list <- list(
  make_option(
    c("-i", "--bedgraph_file"),
    type = "character",
    default = NULL    ,
    metavar = "path"   ,
    help = "beta-values matrix where rows are Sample and columns are length or sample of counts."
  ),
  make_option(
    c("-o", "--outdir"),
    type = "character",
    default = './'    ,
    metavar = "path"   ,
    help = "Output directory."
  )
)


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
bed_data = data.table::fread(opt$bedgraph_file,sep = "\t")
str_bed_name = strsplit(opt$bedgraph_file,"/")
bed_name = str_bed_name[length(str_bed_name)]
out_bed_name = paste(opt$outdir,bed_name,".filter.bed",sep = "")
col_num = dim(bed_data)[2]
new_bed_data = bed_data[rowSums(is.na(bed_data))<col_num*0.6,]
zero_new_bed_data =  new_bed_data[rowSums(new_bed_data[,-c(1,2,3)] == 0,na.rm = TRUE)<20,]
write.table(zero_new_bed_data,file = out_bed_name,quote = FALSE,sep = "\t",col.names = T,row.names = F)
