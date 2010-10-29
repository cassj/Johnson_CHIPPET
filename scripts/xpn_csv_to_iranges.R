#!/usr/bin/Rscript

qw <- function(...) {
  as.character(sys.call()[-1])
}

liftOver<-function(data, chain.file, ucsc.format=T, chr.col="chr", start.col="start",end.col="end"){

  #data should be a matrix or dataframe with cols for chr, start and end
  #TODO: Or a RangedData / IRange object, 
 
  this<-data.frame(chr=as.character(data[,chr.col]),start=as.numeric(data[,start.col]),end=as.numeric(data[,end.col]), stringsAsFactors=F)

  #Normal counting specifies ranges in 1-based, fully closed form.
  #UCSC specifies ranges in  0-based, half open
  if (!ucsc.format){
      this$start<-this$start-1
  }

  #use the chrstartend as an id so we can map back to the original data
  ids = paste(as.character(data[,chr.col]),data[,start.col],data[,end.col], sep=".")
  this <- cbind(this, ids)
  #If we have duplicate positions, remove them for the liftOver 
  this <- unique(this)
  #we also need to chuck out anything that doesn't have positional data.
  no.chr <- which(is.na(this[,1]))
  no.start <- which(is.na(this[,2]))
  no.end <- which(is.na(this[3]))

  inds <- unique(c(no.chr,no.start,no.end))

  if( length(inds)>0 ){ this <- this[-1*inds,] }
  
  
  in.bed <- tempfile()

  #need to watch we don't get scientific notation printed out
  options(scipen=10)
  write.table(this, file=in.bed, sep="\t", row.names=F, col.names=F, quote=F)

  out.bed <- tempfile()
  out.um.bed <- tempfile()
  lo.cmd <- paste("liftOver", in.bed, chain.file, out.bed, out.um.bed)
  system(lo.cmd)

  try(
    new.bed<-read.table(out.bed, sep="\t") 
      ,silent=T)
  try(
    new.um <- read.table(out.um.bed, sep="\t")
      ,silent=T)
  
  #throw away the files
  unlink(c(in.bed, out.bed, out.um.bed))

  if (!exists('new.bed')){stop("No successful mappings")}

  #use the ids as rownames chuck out the column
  #order in the same order as our original data
  #which should stick NAs in for anything that can't be mapped
  rownames(new.bed) <- new.bed[,4]
  new.bed <- new.bed[,-4]
  new.bed <- new.bed[ids,]

  if(!ucsc.format){
   #put the data back to 1-based
   new.bed[,2] <- new.bed[,2]+1
  }

  #replace the new positions in the original dataset
  data[,chr.col] <- new.bed[,1]
  data[,start.col] <- new.bed[,2]
  data[,end.col] <- new.bed[,3]

  #TODO: return some information about the data that won't map
  
  return(data)
  
}


options(stringsAsFactors = FALSE);
args <- commandArgs(trailingOnly=TRUE)

# load datafile from pp_expression_data
limmafile = args[1]
remoat.annot.file = args[2]

limma <- read.csv(limmafile)
limma<-limma[,-1]
colnames(limma)[1]<-"IlluminaID"








######
# Annotation from ReMOAT

# BTW, if we assume that the ~46K probeIDs in the remoat.annot correspond to the entire array
# and that most of those have a 1:1 relationship to a targetID then the "raw" data from KY
# has been seriously filtered. The entire limma table contains only 15060 genes.
# Nothing we can do about it unless we can get hold of the earlier data. Just be wary when
# interpreting results I guess...

# mappings to mm9

# Annotate from the ReMOAT data
# see http://nar.oxfordjournals.org/cgi/content/full/gkp942
remoat.annot <- read.csv(remoat.annot.file, header=T, sep="\t") 

# Throw away anything that isn't perfect or good (the only classes of match likely
# to give a reliable signal according to the ReMOAT paper.
classes <- unique(remoat.annot[,"Quality_score"])
perfect <- grep("Perfect", classes)
good <- grep("Good", classes)
classes <- classes[union(perfect, good)]
keep <- which(remoat.annot[,"Quality_score"] %in% classes)
remoat.annot <- remoat.annot[keep,]

# join illumina data to annotation. This will just ditch any limma data from unmapped probes
# (limma = 15060, limma.annot=13899)
limma.annot = merge(limma, remoat.annot, by.x="IlluminaID", by.y="Target_0")

# The annotation we have is target ID, not probeID.  In most cases the relationship is 1:1
# but there are a few targets which are hit by multiple probes. 
# They seem to be for housekeeping genes (Act1b, Gapdh etc). Presumably controls?
dups<-which(duplicated(limma.annot[,1:7]))
# limma.annot[dups,"IlluminaID"]
rm.inds <- which(limma.annot[,"IlluminaID"] %in% unique(limma.annot[dups,"IlluminaID"]))

# there are only a few, we'll just ditch the duplicates
limma.annot <- limma.annot[-1*rm.inds,]

# We probably don't need *all* the annotation, although we do need to have the genome and
# transcriptome position of the probe.
# Multiple definitions of transcriptome though, so we need results for all of them.
# see http://www.compbio.group.cam.ac.uk/Resources/Annotation/Description_1.0.0.xls
# for colname definitions
limma.annot <- limma.annot[,qw(IlluminaID, logFC, AveExpr, t, P.Value, adj.P.Val, B,Probe_id,Lumi_id,Probe_sequence,Probe_type, Quality_score, Genomic_location,
                               RefSeq_transcripts, Proportion_RefSeq_transcripts,
                               UCSC_transcripts,Proportion_UCSC_transcripts,
                               GenBank_transcripts, Proportion_GenBank_transcripts,
                               Exons,
                               Ensembl_transcripts, Proportion_Ensembl_transcripts,
                               Lumi_transcriptomic_annotation, Lumi_transcriptomic_match )]


#parse the genomic locations. Some span exon junctions and so we'll put
#their 2 genomic locations in as separate entries in the RangedData object...
g.ranges <- limma.annot[,"Genomic_location"]
g.ranges <- sapply(limma.annot[,"Genomic_location"], function(x){strsplit(x,",")})


# split up the ranges where necessary
split.ranges.lengths <- sapply(g.ranges, function(x){length(x)})

# 14 of these have no genomic range - although they hit a transcript.
# None of them look very interesting (low expression, no evidence for \de)
# so I'll just drop them
split.ranges.inds <- rep(1:length(g.ranges), split.ranges.lengths)

# duplicate the rows in the data table where we have split ranges
limma.annot <- limma.annot[split.ranges.inds,]
g.ranges <- unlist(g.ranges)

# parse the g.ranges into "chr","start","end","strand"
g.ranges <- sapply(g.ranges, function(x){strsplit(x,":")})
g.ranges <- do.call("rbind",g.ranges) 
colnames(g.ranges) <- qw(Chr, Start, End, Strand)
top.strand <- which(g.ranges[,"Strand"]=="+")
bottom.strand <- which(g.ranges[,"Strand"]=="-")
g.ranges[top.strand,"Strand"] <- 1
g.ranges[bottom.strand, "Strand"] <- -1

#get rid of the old genomic pos
limma.annot <- limma.annot[,colnames(limma.annot)!=("Genomic_location")]


#create a RangedData object for the region start and end:
library(IRanges)


nms <- paste(g.ranges[,"Chr"], paste(g.ranges[,"Start"],g.ranges[,"End"], sep="-"), g.ranges[,"Strand"], sep=":")

## there are a couple of instances where different (but very similar) probes are mapping
## to the exact same region of the genome. So we can't just use the position as a name:
#wtf <- nms[which(duplicated(nms))]
#wtf <- which(nms %in% wtf)
#limma.annot[wtf,]
#g.ranges[wtf]

nms <- paste(nms," (", limma.annot[,"Probe_id"] ,")", sep="")

#we have to give them names to avoid a bug in ChIPpeakAnnot if we want to use it later
rd.limma <- RangedData(ranges = IRanges(
                         start= as.numeric(g.ranges[,"Start"]),
                         end = as.numeric(g.ranges[,"End"]),
                         names = nms,
                         ),
                       space = g.ranges[,"Chr"],
                       values = cbind(limma.annot, g.ranges[,"Strand"]),
                       universe = "mm9"
                       )

# save the results as RangedData and csv
save(rd.limma, file="limma_rd.R")

rownames(g.ranges) <- rownames(limma.annot) <- nms
write.csv(cbind(g.ranges,limma.annot), file="limma_rd.csv", row.names=F)



