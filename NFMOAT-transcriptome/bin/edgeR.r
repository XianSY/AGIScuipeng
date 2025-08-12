#!/usr/bin/env Rscript
library(edgeR)
library(statmod)
library(ggplot2)
library(optparse)

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

counts = read.table(opt$counts_file,sep = ",",row.names = 1,header = T)
gene_expression = counts
group = read.table(opt$group_file,sep = ",",header = T)
names(group) = c("sample","group")
compare = read.csv(opt$compare_file,header = T)

for(i in seq(1,length(compare[,1]))){
  deg_sample = group[group$group %in% c(compare[i,2],compare[i,3]),1]
  deg_gene_expression = gene_expression[,colnames(gene_expression) %in% deg_sample]
  deg_group = group[group$group %in% c(compare[i,2],compare[i,3]),]
  deg_group$group = factor(deg_group$group)
  deg_gene_expression = deg_gene_expression[,c(deg_group$sample)]
  deglist = DGEList(counts=deg_gene_expression,group = deg_group$group)
  keep = edgeR::filterByExpr(deg_gene_expression,group = group$group)
  deg_gene_expression_keep = deglist[keep,,keep.lib.sizes = FALSE] 
  deglist_norm = calcNormFactors(deg_gene_expression_keep,method = "TMM")
  design = model.matrix(~deg_group$group)
  dge = estimateDisp(deglist_norm,design ,robust = TRUE)
  fit = glmFit(dge,design,robust = TRUE)
  lrt = topTags(glmLRT(fit),n = nrow(deg_gene_expression_keep$counts))
  gene_diff = lrt$table
  gene_diff = as.data.frame(gene_diff)
  gene_diff$type = "normal"
  gene_diff$type[gene_diff$logFC<=-1&gene_diff$FDR<0.05] = "down"
  gene_diff$type[gene_diff$logFC>=1&gene_diff$FDR<0.05] = "up"
  write.csv(gene_diff,file = paste(opt$outdir,compare[i,2],"_vs_",compare[i,3],".csv",sep = ""))
  pic_name = paste(compare[i,1],".png",sep = "")
  pdf_name = paste(compare[i,1],".pdf",sep = "")
  p=ggplot2::ggplot(data = gene_diff,aes(x = logFC,y=-log10(FDR),color=type))+
    geom_point(size = 1)+
    scale_color_manual(values = c('red', 'gray', 'green'), limits = c('up', 'non-significant', 'down')) +  #自定义点的颜色
    labs(x = 'log2 Fold Change', y = '-log10 adjust p-value', title = '', color = '') +  #坐标轴标题
    theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(), #背景色、网格线>、图例等主题修改
          panel.background = element_rect(color = 'black', fill = 'transparent'),
          legend.key = element_rect(fill = 'transparent')) +
    geom_vline(xintercept = c(-1, 1), lty = 3, color = 'black') +  #添加阈值线
    geom_hline(yintercept = -log10(0.05), lty = 3, color = 'black')
  ggplot2::ggsave(pic_name,p,width=7,height=7)
  ggplot2::ggsave(pdf_name,p,width=7,height=7)  
  }
