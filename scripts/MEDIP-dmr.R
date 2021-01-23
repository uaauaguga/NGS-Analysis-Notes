#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(argparse))
parser <- ArgumentParser(description='Differential expression')
parser$add_argument('-i', '--indir', type='character', required=TRUE, help='Dir contain input bam file')
parser$add_argument('-p', '--positive-ids',type='character',required=TRUE, help='Positive samples')
parser$add_argument('-n', '--negative-ids',type='character',required=TRUE, help='Negative samples')
parser$add_argument('-w','--window-size',type='integer',default=50,help='Window size for calculation')
parser$add_argument('--diff-table', type='character',required=TRUE,help='Output difference table')
parser$add_argument('--count-matrix',type='character',required=TRUE,help="Output count matrix")
parser$add_argument('--rpkm-matrix',type='character',required=TRUE,help="Output RPKM matrix")
parser$add_argument('--rms-matrix',type='character',required=TRUE,help="Output rms matrix")
args <- parser$parse_args()

message("Load required packages ...")
suppressPackageStartupMessages(library(MEDIPS))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
message("Done .")


BSgenome <- 'BSgenome.Hsapiens.UCSC.hg38'
uniq <- 0
chr.select <- c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12',
            'chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrM','chrX','chrY')

message('Read positive sample ids')
positive.samples <- read.delim(args$positive_ids,sep="\n",stringsAsFactors=F,header=F)[,1]
message('Read negative sample ids')
negative.samples <- read.delim(args$negative_ids,sep="\n",stringsAsFactors=F,header=F)[,1]
message('Number of positive samples: ', length(positive.samples))
message('Number of negative samples: ', length(negative.samples))

positive.samples.medip <- c()
negative.samples.medip <- c()


message("Load positive samples ...")
for(sample in positive.samples){
  bam.path <- paste0(args$indir,"/",sample,".bam")
  message(sample)
  sample.medips <- MEDIPS.createSet(file = bam.path, BSgenome = BSgenome, uniq = 0, window_size = args$window_size, chr.select = chr.select,paired=TRUE)
  positive.samples.medip  <- c(positive.samples.medip,sample.medips)
  }
message("Done .")
message("Load negative samples ...")
for(sample in negative.samples){
  bam.path <- paste0(args$indir,"/",sample,".bam")
  message(sample)
  sample.medips <- MEDIPS.createSet(file = bam.path, BSgenome = BSgenome, uniq = 0, window_size = args$window_size, chr.select = chr.select,paired=TRUE)
  negative.samples.medip <- c(negative.samples.medip,sample.medips)
}
message("Done .")


message("Perform differential analysis ...")
CS = MEDIPS.couplingVector(pattern = "CG", refObj = sample.medips)
mr.edgeR <- MEDIPS.meth(MSet1 = positive.samples.medip, MSet2 = negative.samples.medip,CSet = CS,diff.method = "edgeR", MeDIP = T, CNV = F, minRowSum = 10)
message("Done .")

message("Post processing ...")
mr.edgeR <- mr.edgeR[!is.na(mr.edgeR$edgeR.p.value),]
mr.edgeR$chr <- trimws(mr.edgeR$chr)
mr.edgeR$start <- trimws(mr.edgeR$start)
mr.edgeR$stop <- trimws(mr.edgeR$stop)
bin.names <- apply(mr.edgeR[,c("chr","start","stop")],1, paste0, collapse = "|",sep="")
names(bin.names) <- NULL
rownames(mr.edgeR) <- bin.names
message("Done .")

cols <- colnames(mr.edgeR)

message("Write count matrix ...")
cols.count <- cols[grep(".bam.counts",cols)]
mat.count <- mr.edgeR[,cols.count]
colnames(mat.count) <- sub(".bam.counts","",cols.count)
write.table(mat.count,file=args$count_matrix,sep='\t',quote=FALSE)
message("Done .")

message("Write RPKM matrix ...")
cols.rpkm <- cols[grep(".bam.rpkm",cols)]
mat.rpkm <- mr.edgeR[,cols.rpkm]
colnames(mat.rpkm) <- sub(".bam.rpkm","",cols.rpkm)
write.table(mat.rpkm,file=args$rpkm_matrix,sep='\t',quote=FALSE)
message("Done .")

message("Write RMS matrix ...")
cols.rms <- cols[grep(".bam.rms",cols)]
mat.rms <- mr.edgeR[,cols.rms]
colnames(mat.rms) <- sub(".bam.rms","",cols.rms)
write.table(mat.rms,file=args$rms_matrix,sep='\t',quote=FALSE)
message("Done .")

message("Write diff table ...")
edger.fields <- c("edgeR.logCPM","edgeR.logFC","edgeR.p.value","edgeR.adj.p.value")
de.table <- mr.edgeR[,edger.fields]
write.table(de.table,file=args$diff_table,sep='\t',quote=FALSE)
message("Done .")


