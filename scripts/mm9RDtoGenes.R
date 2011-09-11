#!/usr/bin/Rscript

options(stringsAsFactors = FALSE);

library(IRanges)

qw <- function(...) {
  as.character(sys.call()[-1])
}

args <- commandArgs(trailingOnly=TRUE)
filename = args[1]
rd <- get(load(filename))


library(ChIPpeakAnno)
library(biomaRt)
ensmart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")


#load mouse transcripts
data(TSS.mouse.NCBIM37)

#The chromosomes names used in teh NCBIM37 dataset are 1..Y, MT and
#various NT_***** non-standard chrs. 

                                        #remove "chr" prefix for  ChIPpeakAnno
names(rd)<-gsub("chr","",names(rd))

#change "M" to "MT" for ChIPpeakAnno
id<-which(names(rd)=="M")
if (length(id)>0){
   names(rd)[id]<-"MT"
}

#chippeakanno discards anything in values, so hang on to them
vals <- as.data.frame(rd)
vals <- vals[,grep("values.", colnames(vals), value=T)]
colnames(vals) <- gsub('values.','',colnames(vals))


# NOTE: TSS.mouse.NCBIM37 is actually *gene* start and end positions,
# not individual transcripts.

#get the most recent annotation data from ensembl
tss <- getAnnotation(ensmart, "TSS")
save(tss, file="tss.RData")


# nearestStart will output the nearest features calculated as peak start - feature start (feature end if feature resides at
#get the nearest start
nearest <- annotatePeakInBatch(rd,
                               AnnotationData=tss,
                               output = "nearestStart" 
                               )



#and later we'll want to know if it overlaps a coding region?
exons <- getAnnotation(ensmart, "Exon") # we should have this file on s3
save(exons, file="exons.RData")

# annotatePeakInBatch is bastard slow, the IRanges findOverlaps function
# works on RangedData and is much faster
overlapping.exon <- findOverlaps(rd, exons, multiple=FALSE)
overlapping.exon<-lapply(overlapping.exon, function(x){
                                 inds <- is.na(x)
                                 x[inds] <- 0
                                 x[!inds] <- 1
                                 return(x)})


overlapping.exon<-do.call(c, overlapping.exon)

#for reasons I don't understand, annotatePeakInBatch doesn't return your
#data in the order you gave it. so, for example rd["1"][1,] is not
#necessarily nearest.tss["1"][1,].

rd.df <- as.data.frame(rd)
ord <- as.character(rd.df$names)

nearest <- as.data.frame(nearest)
rownames(nearest) <- as.character(nearest$peak)
nearest <- nearest[ord,]
nearest <- cbind(nearest,overlapping.exon, vals)



#and that only gives you the ensembl gene ID, so get extra info:                
filters <- c("ensembl_gene_id")
values<-unique(nearest[,"feature"])
attributes <- c("ensembl_gene_id","mgi_symbol", "description")

annot <- getBM(filters=filters, values=values, attributes=attributes, mart=ensmart)

#ditch any that don't have a symbol or a description
no.anno <- intersect(
                     which(annot[,"mgi_symbol"]==""),
                     which(annot[,"description"]==""))
annot <- annot[-1*no.anno,]


# a few have multiple bits of annotation
annot<-cbind(annot, alt.annot="")
dups<-annot[duplicated(annot[,"ensembl_gene_id"]), "ensembl_gene_id"]

#keep the first on and add all the others as alt.annot 
for (d in dups){
  inds <- which(annot[,"ensembl_gene_id"]==d)
  this.alt.annot <- annot[inds[-1], c("mgi_symbol", "description")]
  annot[inds[1],"alt.annot"] <- paste(paste(this.alt.annot[,1], this.alt.annot[,2]), collapse="; ")
}
annot <- annot[!duplicated( annot[,"ensembl_gene_id"] ), ]
rownames(annot) <- annot[,"ensembl_gene_id"]


#add the annotation to your nearest info
nearest <-   cbind(nearest, annot[nearest[,"feature"],   c("mgi_symbol", "description")])


#reorder
colnms <- qw(space,
             start,
             end,
             width,
             names,
             peak,
             strand,
             feature,
             start_position,
             end_position,
             insideFeature,
             distancetoFeature,
             shortestDistance,
             fromOverlappingOrNearest,
             overlapping.exon,
             mgi_symbol,
             description ,
             RE1.Type,
             Left.Motif.Start,
             Left.Motif.End,
             Left.Sequence,
             Right.Motif.End,
             Right.Sequence,
             Right.Direction,
             Spacer
             )
nearest <- nearest[,colnms]

outfile = sub('(.*)\\.R',paste("\\1", '_nearest.csv',sep=""), filename) 
write.csv(nearest, file=outfile, row.names=F)









