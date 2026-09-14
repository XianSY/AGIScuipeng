#!/usr/bin/env Rscript

library(optparse)
library(data.table)

option_list = list(
    make_option(
        c("--input_MatrixeQTL_results","-I"),
        type = "character",
        default = NULL,
        metavar = "path",
        help = "input MatrixeQTL results"
    ),
    make_option(
        c("--snp_info","-P"),
        type = "character",
        default = NULL,
        metavar = "path",
        help = "input SNP postion information"
    ),
    make_option(
        c("--gene_info","-G"),
        type = "character",
        default = NULL,
        metavar = "path",
        help = "input SNP postion information"
    ),
    make_option(
        c("--output","-O"),
        type = "character",
        default = NULL,
        metavar = "path",
        help = "output eQTL results name"
    )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

MatrixeQTL_results = fread(opt$input_MatrixeQTL_results,sep="\t",header=T)
snp_info = fread(opt$snp_info,header=T,sep="\t")
gene_info = fread(opt$gene_info,header=T,sep="\t")
print(colnames(gene_info))
eqtl_merge <- merge( MatrixeQTL_results, snp_info, by.x = "SNP", by.y = "ID_REF", all = FALSE )
eqtl_merge = merge(eqtl_merge,gene_info,by.x = "gene",by.y = "gene_ID",all = FALSE)

setorder(
    eqtl_merge,
    gene,
    CHROM,
    POS
)

eqtl_merge[
    ,
    distance := POS - shift(POS),
    by = .(gene, CHROM)
]

eqtl_merge[
    ,
    cluster := cumsum(
        is.na(distance) | distance > 10000
    ),
    by = .(gene, CHROM)
]

eqtl_merge[
    ,
    cluster_size := .N,
    by = .(gene, CHROM, cluster)
]

eqtl_cluster <- eqtl_merge[
    cluster_size > 3
]

setorderv(
    eqtl_cluster,
    c("gene", "CHROM", "cluster", "p-value"),
    c(1, 1, 1, 1)
)

lead_eqtl <- eqtl_cluster[
    ,
    .SD[1],
    by = .(gene, CHROM, cluster)
]

cluster_position <- eqtl_cluster[
    ,
    .(
        cluster_start = min(POS),
        cluster_end   = max(POS),
        cluster_size  = .N
    ),
    by = .(gene, CHROM, cluster)
]

lead_eqtl <- merge(
    lead_eqtl,
    cluster_position,
    by = c("gene", "CHROM", "cluster"),
    suffixes = c("", "_cluster")
)

lead_eqtl[
    ,
    distance_to_gene := fifelse(
        POS < start,
        start - POS,
        fifelse(
            POS > end,
            POS - end,
            0
        )
    )
]

lead_eqtl[
    ,
    eQTL_type := fifelse(
        distance_to_gene < 1000000,
        "cis",
        "trans"
    )
]

lead_eqtl[
    ,
    c("distance", "cluster_size", "cluster_size_cluster") := NULL
]

setorder(
    lead_eqtl,
    CHROM,
    POS,
    `p-value`
)

fwrite(
    lead_eqtl,
    paste0(opt$output,"eQTL_cluster_significant.txt"),
    sep = "\t",
    quote = FALSE,
    na = "NA"
)
