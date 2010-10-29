#!/usr/bin/Rscript

options(stringsAsFactors = FALSE);
args <- commandArgs(trailingOnly=TRUE)

chain.file <- args[1]

qw <- function(...) {
  as.character(sys.call()[-1])
}


liftOver<-function(data, chain.file, ucsc.format=F, chr.col="chr", start.col="start",end.col="end"){

  #data should be a matrix or dataframe with cols for chr, start and end
  #TODO: Or a RangedData / IRange object, 
 
  this<-data.frame(chr=as.character(data[,chr.col]),start=data[,start.col],end=data[,end.col], stringsAsFactors=F)

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
  
  
  ##all this stuff should be a .C() call but I don't have time to make it work just now.
  in.bed <- tempfile()
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



#read all of the csv files in
as<-read.csv("Altered_Spacer.csv")
cd<-read.csv("Convergent&Divergent.csv")
lhs<-read.csv("Left_Half-Site.csv")
rhs<-read.csv("Right_Half-Site.csv")
canon<-read.csv("Canonical_RE1.csv")
flipped<-read.csv("Flipped_Orientation.csv")
nore1<-read.csv("No_RE1_Motif.csv")

cols<-qw(PET.ID, RE1.Type, Chromosome, 
         Left.Motif.Start, Left.Motif.End, Left.Sequence, Left.Direction, 
         Right.Motif.Start, Right.Motif.End, Right.Sequence, Right.Direction, 
         Spacer)


as<-cbind(as, RE1.Type="Altered_Spacer", stringsAsFactors=F)

#left and right are on the same strand
as<-as[,cols]
data<-as

#this we need to split into 2
converg<-cd[cd[,"Orientation"]=="converge",]
diverg<-cd[cd[,"Orientation"]=="diverge",]

converg<-cbind(converg, RE1.Type="Convergent" )
diverg<-cbind(diverg, RE1.Type="Divergent")

converg<-converg[,cols]
diverg<-diverg[,cols]
data<-rbind(data,converg)
data<-rbind(data, diverg)


colnames(lhs)<- qw(PET.ID, Chromosome, Left.Motif.Start, Left.Motif.End, Left.Sequence, Left.Direction)
lhs<-cbind(lhs, RE1.Type="Left_Half_Site", Right.Motif.Start=NA, Right.Motif.End=NA, Right.Sequence=NA, Right.Direction=NA, Spacer=NA)
data<-rbind(data,lhs)


colnames(rhs)<-qw(PET.ID, Chromosome, Right.Motif.Start, Right.Motif.End, Right.Sequence, Right.Direction)
rhs<-cbind(rhs, RE1.Type="Right_Half_Site", Left.Motif.Start=NA, Left.Motif.End=NA, Left.Sequence=NA, Left.Direction=NA, Spacer=NA)
data<-rbind(data, rhs)


#use LHS for Canonical
colnames(canon)<- qw(PET.ID, Chromosome, Left.Motif.Start, Left.Motif.End, Left.Sequence, Left.Direction)
canon<-cbind(canon, RE1.Type="Canonical", Right.Motif.Start=NA, Right.Motif.End=NA, Right.Sequence=NA, Right.Direction=NA, Spacer=NA )
data<-rbind(data, canon)


flipped<-cbind(flipped, RE1.Type="Flipped_Orientation")
#Left and Right direction always the same
flipped<-flipped[,cols]
data<-rbind(data,flipped)


#fix region start and end before adding sites with no re1
Region.Start<-apply(data,1,
                 function(x){
                    return(
                       min(x[qw(Left.Motif.Start,
                                Left.Motif.End,
                                Right.Motif.Start,
                                Right.Motif.End)],
                       na.rm=T)
		    )
	          }
                )

Region.End<-apply(data,1,
                 function(x){
                    return(
                       max(x[qw(Left.Motif.Start,
                                Left.Motif.End,
                                Right.Motif.Start,
                                Right.Motif.End)],
                       na.rm=T)
		    )
	          }
                )
Region.Start<-as.numeric(Region.Start)
Region.End<-as.numeric(Region.End)

data<-cbind(data, Region.Start=as.numeric(Region.Start), Region.End=as.numeric(Region.End))
cols<-c(cols,qw(Region.Start,Region.End))


nore1<-cbind(nore1, RE1.Type="No_RE1", Left.Motif.Start=NA, Left.Motif.End=NA, Left.Sequence=NA, Left.Direction=NA, Right.Motif.Start=NA, Right.Motif.End=NA, Right.Sequence=NA, Right.Direction=NA, Spacer=NA)
nore1<-nore1[,cols]
data<-rbind(data, nore1)

to.map<-data[,qw(Chromosome, Region.Start, Region.End)]
colnames(to.map) <- c("chr","start","end")
mapped <- liftOver(to.map , chain.file=chain.file)
data[,qw(Region.Start,Region.End)] <- mapped[,qw(start,end)] 

to.map<-data[,qw(Chromosome, Left.Motif.Start, Left.Motif.End)]
colnames(to.map) <- c("chr","start","end")
mapped <- liftOver(to.map , chain.file=chain.file)
data[,qw(Left.Motif.Start, Left.Motif.End)] <- mapped[,qw(start,end)] 

to.map<-data[,qw(Chromosome, Right.Motif.Start, Right.Motif.End)]
colnames(to.map) <- c("chr","start","end")
mapped <- liftOver(to.map , chain.file=chain.file)
data[,qw(Right.Motif.Start,Right.Motif.End)] <- mapped[,qw(start,end)] 



#There's no point in keeping anything that doesn't have a region - we don't know where it maps to:
inds <- union(which(is.na(data$Region.Start)), which(is.na(data$Region.End)))
if(length(inds)>0){
  data <- data[-1*inds,]
}


#create a RangedData object for the region start and end:
library(IRanges)

#we have to give them names to avoid a bug in ChIPpeakAnnot if we want to use it later
rd <- RangedData(ranges = IRanges(
                   start= data$Region.Start,
                   end = data$Region.End,
                  names = as.character(data$PET.ID)
                   ),
                 space = as.character(data$Chromosome),
                 values = data[,
                   qw(RE1.Type,
                      Left.Motif.Start, Left.Motif.End, Left.Sequence, Left.Direction,
                      Right.Motif.Start, Right.Motif.End, Right.Sequence, Right.Direction,
                      Spacer)]

                 )
#And save the result
save(rd, file="RangedData.R")
df <- as.data.frame(rd)
write.csv(df, file="RangedData.csv", row.names=F)
