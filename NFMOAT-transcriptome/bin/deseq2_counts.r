#!/usr/bin/env Rscript
library(ggplot2)
library(optparse)
library(DESeq2)
option_list <- list(
  make_option(
    c("-o", "--outdir"),
    type = "character",
    default = './'    ,
    metavar = "path"   ,
    help = "Output directory."
  ),
   make_option(
    c("-C", "--counts_file"),
    type = "character",
    default = NULL    ,
    metavar = "path"   ,
    help = "Counts matrix where rows are gene and columns are sample."
  ),
  make_option(
    c("-g", "--group_file"),
    type = "character",
    default = NULL    ,
    metavar = "path"   ,
    help = "Group matrix where rows are Sample and columns are group ."
  ),
  make_option(
    c("-c", "--compare_file"),
    type = "character",
    default = NULL    ,
    metavar = "path"   ,
    help = "Compare matrix where rows are Sample and columns are comparation."
  )
)


opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

group = read.csv(opt$group_file,header = T)
compare = read.csv(opt$compare_file,header = T)    

gene_exs = read.csv(opt$counts_file,header = T,sep=",",row.names = 1)
gene_expression = gene_exs

DEG_count = data.frame(
  DEG_set = character(),
  DEG_number = integer(),
  up_regulated = integer(),
  down_regulated = integer()
)

for(i in seq(1,length(compare[,1]))){
  deg_sample = group[group$group %in% c(compare[i,2],compare[i,3]),1]
  deg_gene_expression = gene_expression[,colnames(gene_expression) %in% deg_sample]
  deg_group = group[group$group %in% c(compare[i,2],compare[i,3]),]
  dds = DESeq2::DESeqDataSetFromMatrix(countData = deg_gene_expression,colData = deg_group,design = ~ group)  
  deg = DESeq2::DESeq(dds)
  res = DESeq2::results(deg)
  results = as.data.frame(res)
  results = results[!is.na(results$log2FoldChange),]
  results$gene_type = NA
  results$gene_type[results$log2FoldChange>=1 & results$pvalue<0.05] = "up"
  results$gene_type[results$log2FoldChange<=-1 & results$pvalue<0.05] = "down"
  write.csv(results,file = paste(opt$outdir,compare[i,2],"_vs_",compare[i,3],".csv",sep = ""),quote = F)
  test_DEG = data.frame(
  DEG_set = compare[i, 1],
  DEG_number = sum(!is.na(results$gene_type)),
  up_regulated = sum(results$gene_type == 'up', na.rm = TRUE),
  down_regulated = sum(results$gene_type == "down", na.rm = TRUE)
  )
  DEG_count = rbind(DEG_count, test_DEG)

  png_name = paste(compare[i,1],".png",sep = "")
  pdf_name = paste(compare[i,1],".pdf",sep = "")
  pl = ggplot2::ggplot(data = results, aes(x = log2FoldChange, y = -log10(pvalue))) + 
    geom_point(aes(color = gene_type),alpha = 1,shape = 19,size = 1.5) +
    theme_bw() + theme(panel.grid.minor = element_blank(),
                       panel.grid.major = element_blank()) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(log2(0.5), log2(2)), linetype = "dashed")
  ggplot2::ggsave(png_name,pl,width=7,height=7) 
  ggplot2::ggsave(pdf_name,pl,width=7,height=7)
 }
write.csv(DEG_count,file = paste(opt$outdir,"DEG_counts",".csv",sep = ""),quote = F,row.names = F)

