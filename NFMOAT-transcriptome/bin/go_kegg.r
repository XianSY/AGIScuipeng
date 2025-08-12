#!/usr/bin/env Rscript
library(ggplot2)
library(optparse)
library(clusterProfiler)

option_list = list(
    make_option(
        c("-i","--input_gene"),
        type = "character",
        default = NULL,
        metavar = "path",
        help = "please input gene list"
    ),
    make_option(
        c("-S","--species"),
        type = "character",
        default = NULL,
        metavar = "STRING",
        help = "please input species name"
    )
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

enrichGene = read.table(opt$input_gene,header=F)
names(enrichGene) = "gene"

GO = function(gene,GOannotation){
  GOannotation = split(GO_backgroud,with(GO_backgroud,ontology))
  term2gene = GOannotation$BP[,c("goid","gene_id")]
  term2name = GOannotation$BP[,c("goid","desc")]
  ego = clusterProfiler::enricher(
    gene = gene,
    TERM2GENE = term2gene,
    TERM2NAME = term2name
  )
  BP = ego@result
  BP$ontology = "BP"
  
  term2gene = GOannotation$MF[,c("goid","gene_id")]
  term2name = GOannotation$MF[,c("goid","desc")]
  ego = clusterProfiler::enricher(
    gene = gene,
    TERM2GENE = term2gene,
    TERM2NAME = term2name
  )
  MF= ego@result
  MF$ontology = "MF"
  
  term2gene = GOannotation$CC[,c("goid","gene_id")]
  term2name = GOannotation$CC[,c("goid","desc")]
  ego = clusterProfiler::enricher(
    gene = gene,
    TERM2GENE = term2gene,
    TERM2NAME = term2name
  )
  CC= ego@result
  CC$ontology = "CC"
  
  GO_results = rbind(BP,MF,CC)
  return(GO_results)
}

KEGG = function(gene,KEGG_backgroud){
  term2gene = KEGG_backgroud[,c("ID","geneID")]
  term2name = KEGG_backgroud[,c("ID","Description")]
  kegg <- enricher( gene = gene,
                    TERM2GENE = term2gene, 
                    TERM2NAME = term2name,)
  kegg_results = kegg@result
  return(kegg_results)
}

if( opt$species == "Oryza"){
    GO_backgroud  = read.table("/opt/data/rice_GO.csv",header=T,sep = ";")
    KEGG_backgroud = read.table("/opt/data/rice_KEGG.csv",header=T, sep=";")
    GO_results = GO(enrichGene$gene,GO_backgroud)
    KEGG_results = KEGG(enrichGene$gene,KEGG_backgroud)

    write.csv(GO_results,file = paste(opt$species,"GO_results.csv",sep = "_"),quote = F,col.names = T,row.names = F)
    write.csv(KEGG_results,file = paste(opt$species,"KEGG_results.csv",sep = "_"),quote = F,col.names = T,row.names = F)
}else if(opt$species == "Soybean"){
    GO_backgroud  = read.table("/opt/data/soybean_GO.csv",header=T,sep = ";")
    KEGG_backgroud = read.table("/opt/data/soybean_KEGG.csv",header=T, sep="=")
    GO_results = GO(enrichGene$gene,GO_backgroud)

    term2gene = KEGG_backgroud[,c("KO","geneid")]
    term2name = KEGG_backgroud[,c("KO","Level_2")]
    kegg <- enricher( gene = enrichGene$gene,
                    TERM2GENE = term2gene, 
                    TERM2NAME = term2name,)
    KEGG_results = kegg@result

    write.csv(GO_results,file = paste(opt$species,"GO_results.csv",sep = "_"),quote = F,col.names = T,row.names = F)
    write.csv(KEGG_results,file = paste(opt$species,"KEGG_results.csv",sep = "_"),quote = F,col.names = T,row.names = F)
}else if(opt$species == "Zea"){
    GO_backgroud  = read.table("/opt/data/zea_GO.csv",header=T,sep = ";")
    KEGG_backgroud = read.table("/opt/data/zea_kegg.csv",header=T, sep="=")
    GO_results = GO(enrichGene$gene,GO_backgroud)
    KEGG_results = KEGG(enrichGene$gene,KEGG_backgroud)
    write.csv(GO_results,file = paste(opt$species,"GO_results.csv",sep = "_"),quote = F,col.names = T,row.names = F)
    write.csv(KEGG_results,file = paste(opt$species,"KEGG_results.csv",sep = "_"),quote = F,col.names = T,row.names = F)
}



