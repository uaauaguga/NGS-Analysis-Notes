#!/usr/bin/env Rscript
library(GenomicFeatures)
gtf.path <- "Arabidopsis_thaliana.TAIR10.34.gtf"

message("Load gff file ...")
txdb <- makeTxDbFromGFF(file=gtf.path,format="gtf")


message("Get 5pUTR ...")
tx.utr5p <- fiveUTRsByTranscript(txdb,use.names=TRUE)


message("Get CDS ...")
tx.cds <- cdsBy(txdb,by="tx",use.names=TRUE)

message("Get exons ...")
tx.exons <- exonsBy(txdb, by="tx",use.names=TRUE)


message("Get location ...")

utr5p.lengths <- sum(width(tx.utr5p))
cds.lengths <- sum(width(tx.cds))
tx.lengths <- sum(width(tx.exons))

tx.no5putr <- setdiff(names(cds.lengths),names(utr5p.lengths))
utr5p.add <- rep(0,length(tx.no5putr))
names(utr5p.add) <- tx.no5putr
utr5p.lengths <- c(utr5p.lengths,utr5p.add)
utr5p.lengths <- utr5p.lengths[names(cds.lengths)]
tx.lengths <- tx.lengths[names(cds.lengths)]

tx.coordinates <- data.frame(utr5p.lengths,cds.lengths,tx.lengths)
tx.coordinates[["utr5p-start"]] <- 0
tx.coordinates[["utr5p-end"]] <- tx.coordinates[["utr5p.lengths"]]
tx.coordinates[["cds-end"]] <- tx.coordinates[["utr5p.lengths"]] + tx.coordinates[["cds.lengths"]]
tx.coordinates[["utr3p-end"]] <- tx.coordinates[["tx.lengths"]]

utr5p <- tx.coordinates[,c("utr5p-start","utr5p-end")]
utr5p <- utr5p[utr5p[["utr5p-end"]]>0,]

utr3p <- tx.coordinates[,c("cds-end","utr3p-end")]
print(utr3p)
utr3p <- utr3p[utr3p[["cds-end"]]<utr3p[["utr3p-end"]],]

message("Save location ...")

write.table(utr5p,file="tair10.tx.utr5p.bed",sep="\t",quote = F, col.names = F)
write.table(tx.coordinates[,c("utr5p-end","cds-end")],file="tair10.tx.cds.bed",sep="\t",quote = F, col.names = F)
write.table(utr3p,file="tair10.tx.utr3p.bed",sep="\t",quote = F, col.names = F)


message("All done .")