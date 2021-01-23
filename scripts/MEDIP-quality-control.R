#!/usr/bin/env Rscript
message("Load required packages ...")
suppressPackageStartupMessages(library(MEDIPS))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressPackageStartupMessages(library(argparse))
message("Done .")

parser <- ArgumentParser(description='Differential expression')
parser$add_argument('-i', '--input', type='character', required=TRUE,help='Input bam file')
parser$add_argument('-m', '--method', type='character',default='enrichment',choices=c('enrichment','saturation'),help='Operation to perform')
parser$add_argument('-w','--window-size',type='integer',default=50,help='Window size for calculation')
parser$add_argument('-o', '--output', type='character',required=TRUE,help='Output')
args <- parser$parse_args()


BSgenome <- 'BSgenome.Hsapiens.UCSC.hg38'
uniq <- 0
chr.select <- c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12',
            'chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrM','chrX','chrY')
if(args$method=='enrichment'){
  er = MEDIPS.CpGenrich(file = args$input, BSgenome = BSgenome, chr.select = chr.select,uniq = uniq,paired=TRUE)
  #write.table(er, file=args$output, sep='\t',row.names = FALSE,quote=FALSE)
  er <- t(data.frame(er))
  write.table(er,file=args$output,sep="\t",col.names=FALSE,quote=FALSE)
}else if(args$method=='saturation'){
  sr = MEDIPS.saturation(file = args$input, BSgenome = BSgenome, chr.select = chr.select, uniq = uniq, window_size = args$window_size,paired=TRUE)
  saturation.table <- sr$distinctSets
  colnames(saturation.table) <- c('subset','correlation')
  saturation.table <- transform(saturation.table,data=rep("observed",nrow(saturation.table)))
  estimate.table <- sr$estimation
  colnames(estimate.table) <- c('subset','correlation')
  estimate.table <- transform(estimate.table,data=rep("estimated",nrow(estimate.table)))
  saturation.table <-rbind(saturation.table,estimate.table)
  write.table(saturation.table, file=args$output, sep='\t',row.names = FALSE,quote=FALSE)
}
