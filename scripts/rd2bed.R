#!/usr/bin/Rscript

options(stringsAsFactors = FALSE);

library(IRanges)

qw <- function(...) {
  as.character(sys.call()[-1])
}

args <- commandArgs(trailingOnly=TRUE)
filename = args[1]
filename <- "/mnt/data/NS5/PET/RangedData.R"
outfile <- sub('RangedData.R','peaks.bed', filename)

rd <- get(load(filename))

rd <- as.data.frame(rd)
bed <- rd[,qw(space, start,end)]
bed<- bed[order(rd[,"space"],rd[,"start"]),]

write.table(bed, file=outfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
