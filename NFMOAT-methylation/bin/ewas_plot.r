#!/usr/bin/env Rscript

library(CMplot)
library(dplyr)
library(tidyr)
library(optparse)

#读取数据
option_list = list(
  make_option(
    c("-i","--input_ewas_results"),
    type = "character",
    default = NULL,
    metavar = "path",
    help = "phenotype csv file "
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


#读取ewas结果数据
ewas_res = read.table(opt$input_ewas_results, header = TRUE, sep = "")
ewas_res_plot = ewas_res %>% 
  dplyr::mutate(ewas_site = Probe) %>%
  tidyr::separate(Probe, into = c("chr", "BPs")) %>% 
  select(c("ewas_site", "chr", "BPs", "p"))

str_bed_name = strsplit(opt$input_ewas_results,"/")
bed_name = str_bed_name[length(str_bed_name)]
out_bed_name = paste0(bed_name,".plot")

#保存画图结果数据
write.table(ewas_res_plot,file = out_bed_name,sep = "\t",quote = FALSE,row.names = FALSE)

#QQ plot
CMplot(
  ewas_res_plot,
  plot.type = "q",
  threshold = 0.05,
  file.name = out_bed_name,
)

#曼哈顿图
CMplot(
  ewas_res_plot,
  plot.type = "m",
  threshold = c(0.01, 0.05) / nrow(ewas_res_plot),
  threshold.col = c('grey', 'black'),
  threshold.lty = c(1, 2),
  threshold.lwd = c(1, 1),
  amplify = T,
  file.name = out_bed_name,
  signal.cex = c(1, 1),
  signal.pch = c(20, 20),
  signal.col = c("red", "orange")
)

CMplot(
  ewas_res_plot,
  plot.type = "m",
  threshold = c(0.01, 0.05) / nrow(ewas_res_plot),
  threshold.col = c('grey', 'black'),
  threshold.lty = c(1, 2),
  threshold.lwd = c(1, 1),
  amplify = T,
  file="pdf",
  file.name = out_bed_name,
  signal.cex = c(1, 1),
  signal.pch = c(20, 20),
  signal.col = c("red", "orange")
)

#染色体密度图
CMplot(
  ewas_res_plot,
  plot.type = "d",
  bin.size = 1e6,
  file.name = out_bed_name,
  col = c("darkgreen", "yellow", "red"),
)

CMplot(
  ewas_res_plot,
  plot.type = "d",
  bin.size = 1e6,
  file = "pdf",
  file.name = out_bed_name,
  col = c("darkgreen", "yellow", "red"),
)


#环形曼哈顿图
CMplot(
  ewas_res_plot,
  plot.type = "c",
  r = 0.4,
  col = c("grey30", "grey60"),
  threshold = c(1e-6, 1e-4),
  cir.chr.h = 1.5,
  amplify = TRUE,
  threshold.lty = c(1, 2),
  threshold.col = c("red", "blue"),
  signal.line = 1,
  signal.col = c("red", "green"),
  chr.den.col = c("darkgreen", "yellow", "red"),
  bin.size = 1e6,
  outward = FALSE,
  dpi = 300,
  file.output = TRUE,
  verbose = TRUE,
  file.name = out_bed_name,
)

#环形曼哈顿图
CMplot(
  ewas_res_plot,
  plot.type = "c",
  r = 0.4,
  col = c("grey30", "grey60"),
  threshold = c(1e-6, 1e-4),
  cir.chr.h = 1.5,
  amplify = TRUE,
  threshold.lty = c(1, 2),
  threshold.col = c("red", "blue"),
  signal.line = 1,
  signal.col = c("red", "green"),
  chr.den.col = c("darkgreen", "yellow", "red"),
  bin.size = 1e6,
  outward = FALSE,
  dpi = 300,
  file.output = TRUE,
  verbose = TRUE,
  file="pdf",
  file.name = out_bed_name,
)

