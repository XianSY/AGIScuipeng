#!/opt/conda/envs/emmax-anno/bin/Rscript --verbose
library(ChIPseeker)
library(GenomicFeatures)

args <- commandArgs(TRUE)

bedgraph = args[1]
gtf_file = args[2]

names = strsplit(bedgraph,"/")

gtf = GenomicFeatures::makeTxDbFromGFF("Oryza_sativa.gff")
T1 = ChIPseeker::readPeakFile(
  bedgraph
)

peakAnno_T1 = ChIPseeker::annotatePeak(T1,tssRegion = c(-3000,3000),TxDb = gtf)
peakAnno_df = as.data.frame(peakAnno_T1)

write.csv(peakAnno_df,file = paste(names,".csv",sep = ""))

png("pie.png",width = 600,height = 800)
ChIPseeker::plotAnnoPie(peakAnno_T1)
dev.off()

png("TSS.png",width = 600,height = 800)
ChIPseeker::plotAnnoPie(peakAnno_T1)
dev.off()
