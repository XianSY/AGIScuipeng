#!/usr/bin/env Rscript

library(optparse)

option_list = list(
  make_option(
    c("--beta_matrix","-E"),
    default = NULL,
    type = "character",
    metavar = "path",
    help = "expression files"
  )
)


opt_parser = OptionParser(option_list = option_list)
opt        = parse_args(opt_parser)
bedgraph = read.table(opt$beta_matrix,header = T)

bedgraph$smpID=paste(bedgraph$chrom,bedgraph$start,sep = ":")
new_bedgraph = bedgraph[,c(length(bedgraph),seq(1,length(bedgraph)-1))]

smp_bed = new_bedgraph[,c("smpID","chrom","start","end")]
smp_beta = new_bedgraph[,c(1,seq(5,length(new_bedgraph)))]

write.table(smp_bed,"smp_pos.txt",quote = F,col.names = T,row.names = F,sep = "\t")
write.table(smp_beta,"smp_beta.txt",quote = F,col.names = T,row.names = F,sep = "\t")
