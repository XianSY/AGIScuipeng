#!/usr/bin/env Rscript



library(optparse)

option_list = list(
  make_option(
    c("--bed","-B"),
    type = "character",
    default = NULL,
    metavar = "path",
    help = "gene bed informations"
  ),
  make_option(
    c("--wgt_pos","-W"),
    type = "character",
    default = NULL,
    metavar = "path",
    help = "weights rdata statics file "
  )
)

opt_parse = OptionParser(option_list = option_list)
opt       = parse_args(opt_parse)

if(is.null(opt$bed)){
  print_help(opt_parser)
  stop("Please provide a gene bed informations file .",call. = FALSE)
}else if(is.null(opt$wgt_pos)){
  print_help(opt_parser)
  stop("Please provide a weights rdata statics file .",call. = FALSE)
}

gene_info = read.table(opt$bed,header = F)
names(gene_info) = c("ID","CHR","P0","P1")

gene_wgt_info = read.table(opt$wgt_pos,header = F)
names(gene_wgt_info) = c("WGT","ID")

wegiht_pos = merge(gene_wgt_info,gene_info,by = "ID")
wegiht_pos = wegiht_pos[,c("WGT","ID","CHR","P0","P1")]
write.table(wegiht_pos,file = "gene_wgt.pos",quote = F,row.names = F,col.names = T,sep = "\t")
