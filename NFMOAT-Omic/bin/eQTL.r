#!/usr/bin/env Rscript

library(optparse)
library(MatrixEQTL)
library(data.table)
#library(parallel)
options(mc.cores=16)




option_list = list(
  make_option(
    c("--gene_info","-B"),
    default = NULL,
    metavar = "path",
    type    = "character",
    help    = "gene information of bed"
  ),
  make_option(
    c("--snp_pos","-S"),
    default = NULL,
    metavar = "path",
    type    = "character",
    help    = "snp position",
  ),
  make_option(
    c("--snp_vcf_num"),
    default = NULL,
    type    = "character",
    metavar = "path",
    help    = "vcf .raw information" 
  ),
  make_option(
    c("--gene_expression"),
    default = NULL,
    type    = "character",
    metavar = "path",
    help    = "gene expression of TPM(FPKM RPKM)" 
  ),
  make_option(
    c("--output_name"),
    default = NULL,
    type    = "character",
    metavar = "path",
    help    = "please input output_name"
)
)
opt_parse = OptionParser(option_list = option_list)
opt       = parse_args(opt_parse)


gene_pos = data.table::fread(opt$gene_info,sep = "\t",header=TRUE)
snp_pos = data.table::fread(opt$snp_pos,sep = "\t",header=TRUE)
snps = SlicedData$new()
snps$fileDelimiter = "\t"       #指定SNP文件的分隔符
snps$fileOmitCharacters = "NA" #定义缺失值
snps$fileSkipRows = 1        #跳过第一行（适用于第一行是列名的情况）
snps$fileSkipColumns = 1     #跳过第一列（适用于第一列是SNP ID的情况）
snps$fileSliceSize = 20000      #每次读取2000条数据
snps$LoadFile(opt$snp_vcf_num)



gene = SlicedData$new()
gene$fileDelimiter = "\t" 
gene$fileOmitCharacters = "NA"
gene$fileSkipRows = 1
gene$fileSkipColumns = 1
gene$fileSliceSize = 2000
gene$LoadFile(opt$gene_expression)


output_file_name_cis = paste0(opt$output_name,"_","cis_eqtl.txt") #顺式
output_file_name_tra = paste0(opt$output_name,"_","trans_eqtl.txt") #反式 # 只有在高于阈值重要的gene-SNP关联对才会被保存
pvOutputThreshold_cis = 0.05/(as.numeric(nrow(snp_pos)))
pvOutputThreshold_tra = 0.05/(as.numeric(nrow(snp_pos)))

useModel = modelLINEAR # modelANOVA, modelLINEAR, or modelLINEAR_CROSS 

my_eQTL = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  useModel = useModel,
  verbose = TRUE,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra, 
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis, 
  snpspos = snp_pos,
  genepos = gene_pos,
  cisDist = 1e6,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE
)


out = my_eQTL$all$eqtls
write.table(out,file = "eqtlout.txt")
pdf("p_value_qqplots.pdf")
plot(my_eQTL)
dev.off()

q()

