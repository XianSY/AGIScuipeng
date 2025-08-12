#!/usr/bin/env Rscript

library(optparse)


option_list = list(
  make_option(
    c("--gwas_resluts","-R"),
    type = "character",
    default = NULL,
    metavar = "path",
    help = "emmax software to GWAS results",
  ),
  make_option(
    c("--allele_info","-A"),
    type = "character",
    default = NULL,
    metavar = "path",
    help = "allele gene informations (plink to gain)"
  )
)

opt_parse = OptionParser(option_list = option_list)
opt       = parse_args(opt_parse)


gwas_resluts = read.table(opt$gwas_resluts,header = F)
names(gwas_resluts)[1] = "SNP"


allele_info = read.table(opt$allele_info,header = T)

LD_sore = merge(gwas_resluts,allele_info,by = "SNP")
LD_sore$Z = LD_sore$V2/LD_sore$V3
LD_sore = LD_sore[,c("SNP","A1","A2","Z")]
write.table(LD_sore,file = "LD.sumstats",quote = F,col.names = T,row.names = F)
