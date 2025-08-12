#!/usr/bin/env Rscript
library(optparse)
library(data.table)

option_list = list(
  make_option(
    c("--vcf_num","-E"),
    default = NULL,
    type = "character",
    metavar = "path",
    help = "expression files"
  )
)


opt_parser = OptionParser(option_list = option_list)
opt        = parse_args(opt_parser)

if(is.null(opt$vcf_num)){
  print_help(opt_parser)
  stop("Please provide a vcf_number file .",call. = FALSE)
}

convert_genotype <- function(gt) {
  sapply(gt, function(gt_field) {
    # 提取 GT 字段
    gt <- strsplit(gt_field, ":")[[1]][1]
    
    # 替换基因型格式
    gt <- gsub("0/0|0\\|0", "0", gt)
    gt <- gsub("0/1|0\\|1|1\\|0", "1", gt)
    gt <- gsub("1/1|1\\|1", "2", gt)
    gt <- gsub("\\./\\.|\\.\\|\\.", "NA", gt)
    
    # 转换为数值类型
    return(as.numeric(gt))
  })
}

vcf_lines <- readLines(opt$vcf_num)
header_line <- vcf_lines[grep("^#CHROM", vcf_lines)]
data_lines <- vcf_lines[!grepl("^#", vcf_lines)]
# 处理头部
header <- unlist(strsplit(header_line, "\t"))
# 处理数据行
data_matrix <- do.call(rbind, strsplit(data_lines, "\t"))
colnames(data_matrix) <- header
# 创建数据框
geno_data <- as.data.frame(data_matrix, stringsAsFactors = FALSE)
# 提取基因型列（从 FORMAT 列开始）
format_col_index <- which(header == "FORMAT")
geno_cols <- (format_col_index + 1):ncol(geno_data)
# 转换基因型
genotypes_list <- list()
for(i in geno_cols) { 
  genotypes_list[[i-format_col_index]] <- convert_genotype(geno_data[,i])}
# 合并成矩阵
genotypes <- do.call(cbind, genotypes_list)
rownames(genotypes) <- paste0(geno_data[,1], ":", geno_data[,2])  # CHROM_POS
colnames(genotypes) <- header[geno_cols]
write.table(genotypes,file = "popu_snp.chromesome_t.txt",col.names = T,row.names = T,quote = F,sep = "\t")
# 创建SNP位置文件
snps_location <- data.frame(  
  snpid = rownames(genotypes),  
  chr   = geno_data[,1],  # CHROM column  
  start = as.numeric(geno_data[,2])-1,
  end   = as.numeric(geno_data[,2]),  # POS column
  stringsAsFactors = FALSE)
write.table(snps_location,file = "snps_pos.txt",col.names = T,row.names = F,quote = F,sep = "\t")
