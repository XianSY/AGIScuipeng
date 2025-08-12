#!/usr/bin/env Rscript

library(optparse)
library(ggplot2)
library(data.table)

option_list = list(
  make_option(
    c("-o","--outdir"),
    type = "character",
    default = "./",
    metavar = "path",
    help = "Output directory"
  ),
  make_option(
    c("--bedgraph"),
    type = "character",
    default = NULL,
    metavar = "path",
    help = "Input directory of bedgraph"
  ),
  make_option(
    c("--cx"),
    type = "character",
    default = NULL,
    metavar = "value",
    help = "Input context of bedgraph"
  )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

cx = opt$cx

if(cx=="cpg"){
  cx_col = "#E41A1C"
} else if (cx == "chg"){
  cx_col = "#377EB8"
}else if (cx=="chh"){
  cx_col = "#4DAF4A"
}

bedgraph = data.table::fread(input = opt$bedgraph,header = TRUE,sep = "\t")
one_three = bedgraph[,c(1,2,3)]
for(i in names(bedgraph)[seq(4,length(bedgraph))]){
  sample = i  
  sub_bedgraph = cbind(bedgraph[,1:3],bedgraph[[i]])
    names(sub_bedgraph) = c(names(bedgraph)[seq(1,3)],"level")
    sub_bedgraph = sub_bedgraph[sub_bedgraph$level != 0, ]
    p = ggplot2::ggplot(data = sub_bedgraph, aes(x = level)) + theme_minimal() + geom_histogram(
      aes(y = after_stat(density)),
      alpha = 0.9,
      col = "black",
      fill = cx_col,
      binwidth = 11,
      boundary = 10
    ) + theme_bw()
  hist_name = paste(sample,".png",sep = "")
  print(hist_name)
  ggplot2::ggsave(hist_name,p,width=7,height=7)
}
