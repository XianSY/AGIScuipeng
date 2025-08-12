#!/usr/bin/envs Rscript

library(optparse)

#读取数据
option_list = list(
  make_option(
    c("-i","--dmr_chg_results"),
    type = "character",
    default = NULL,
    metavar = "path",
    help = "chg dmr results"
  ),
  make_option(
    c("-i","--dmr_cpg_results"),
    type = "character",
    default = NULL,
    metavar = "path",
    help = "cpg dmr results"
  ),
  make_option(
    c("-i","--dmr_chh_results"),
    type = "character",
    default = NULL,
    metavar = "path",
    help = "chh dmr results"
  ),
  make_option(
    c("-o","--outdir"),
    type = "character",
    default = "./",
    metavar = "path",
    help = "Output phenotype file"
  )
)

#解析参数
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

cpg_header = c(
  "chr",
  "start",
  "stop",
  "q-value",
  "mean methylation difference",
  "#CpGs",
  "p (MWU)",
  "p (2D KS)",
  "mean g1",
  "mean g2"
)

chg_header = c(
  "chr",
  "start",
  "stop",
  "q-value",
  "mean methylation difference",
  "#CHGs",
  "p (MWU)",
  "p (2D KS)",
  "mean g1",
  "mean g2"
)
chh_header = c(
  "chr",
  "start",
  "stop",
  "q-value",
  "mean methylation difference",
  "#CHHs",
  "p (MWU)",
  "p (2D KS)",
  "mean g1",
  "mean g2"
)




#CHG_dmr
chg_dmr_results = read.table(opt$dmr_chg_results,header = F)
chg_dmr_results_name = strsplit(opt$dmr_chg_results,"/")
chg_dmr_name = chg_dmr_results_name[length(chg_dmr_results_name)]
colnames(chg_dmr_results) = chg_header

vs_names_strsplit = strsplit(chg_dmr_results_name = strsplit,"_")


chg_hyper = chg_dmr_results[chg_dmr_results$`mean methylation difference`>0 & chg_dmr_results$`p (2D KS)`<0.05,]
chg_hypo = chg_dmr_results[chg_dmr_results$`mean methylation difference`< 0  & chg_dmr_results$`p (2D KS)`<0.05,]
write.table(chg_hyper,file = paste0(chg_dmr_name,"chg_hyper.txt"),quote = F,col.names = T,row.names = F)
write.table(chg_hypo,file = paste0(chg_dmr_name,"chg_hypo.txt"),quote = F,col.names = T,row.names = F)


#CpG_dmr
cpg_dmr_results = read.table(opt$dmr_cpg_results,header = F)
colnames(cpg_dmr_results) = cpg_header
cpg_dmr_results_name = strsplit(opt$dmr_cpg_results,"/")
cpg_dmr_name = cpg_dmr_results_name[length(cpg_dmr_results_name)]


cpg_hyper = cpg_dmr_results[cpg_dmr_results$`mean methylation difference`>0 & cpg_dmr_results$`p (2D KS)`<0.05,]
cpg_hypo = cpg_dmr_results[cpg_dmr_results$`mean methylation difference`< 0  & cpg_dmr_results$`p (2D KS)`<0.05,]
write.table(cpg_hyper,file = paste0(cpg_dmr_name,"cpg_hyper.txt"),quote = F,col.names = T,row.names = F)
write.table(cpg_hypo,file = paste0(cpg_dmr_name,"cpg_hypo.txt"),quote = F,col.names = T,row.names = F)


#CHH_dmr
chh_dmr_results = read.table(opt$dmr_chh_results,header = F)
colnames(chh_dmr_results) = chh_header
chh_dmr_results_name = strsplit(opt$dmr_chh_results,"/")
chh_dmr_name = chh_dmr_results_name[length(chh_dmr_results_name)]

chh_hyper = chh_dmr_results[chh_dmr_results$`mean methylation difference`>0 & chh_dmr_results$`p (2D KS)`<0.05,]
chh_hypo = chh_dmr_results[chh_dmr_results$`mean methylation difference`< 0  & chh_dmr_results$`p (2D KS)`<0.05,]
write.table(cpg_hyper,file = paste0(chh_dmr_name,"cpg_hyper.txt"),quote = F,col.names = T,row.names = F)
write.table(cpg_hypo,file =  paste0(chh_dmr_name,"cpg_hypo.txt"),quote = F,col.names = T,row.names = F)



dmr_results = data.frame(hyper = c(dim(chg_hyper)[1], dim(cpg_hyper)[1], dim(chh_hyper)[1]), hypo = c(dim(cpg_hypo)[1], dim(chg_hypo)[1], dim(chh_hypo)[1]), row.names = c("chg", "cpg", "chh"))

write.table(dmr_results,file = paste0(vs_names_strsplit[4],vs_names_strsplit[6],"dmr_results.txt",quote = F,col.names = T,row.names = F)

