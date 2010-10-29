#!/usr/bin/Rscript

# because I can't be bothered to type all the quotes
qw <- function(...) {
  as.character(sys.call()[-1])
}

options(stringsAsFactors = FALSE);
options(scipen=10)

args <- commandArgs(trailingOnly=TRUE)
filename = args[1]
outfile = args[2]
data <- read.csv(filename)

# This data seems to have been pre-filtered and bg-corrected
# It isn't clear how this has been done, and there doesn't seem to
# be any way of getting hold of the original data.

#grab expression value cols
cols <- qw(dn, dn, dn, dn, ev, ev, ev, ev)
rownames(data)<-data[,1]
n<-length(cols)
data<-data[,2:(n+1)]
ev<-which(cols=="ev")
dn<-which(cols=="dn")


# set values <10 to 10
# This is a kludge and could possibly have horrible repercussions later on when
# we're determining differential expression as there will be some cases in which there
# is 0 variation. Not really sure what else to do. Limma fudges the variance anyway, so
# hopefully it should be reasonably robust to this sort of thing, but even so...
# Should be fine for just getting an idea of the things that are *really* changing, but
# be careful using this dataset for anything more complicated without maybe getting hold of
# the raw data.
fix10<-function(x){
  x[x<10]<-10
  return(x)
}

data<-apply(data,2,fix10)

#quantile normalise. The normaliseIllumina method seems to cope with 
#an expressionset, which save the hassle of trying to build a valid
#BSData object
library(affy)
library(beadarray)

es <- new("ExpressionSet", exprs = data)
E = normaliseIllumina(es, method="quantile", transform="log2")
data <- exprs(E)

library(limma)

design<-matrix(0,nrow=(ncol(data)), ncol=2)
colnames(design)<-c("ev","dn")
design[ev,"ev"]<-1
design[dn,"dn"]<-1


fit<-lmFit(data, design)
cont.matrix<-makeContrasts(dnvsev=dn-ev, levels=design)
fit<-contrasts.fit(fit, cont.matrix)
ebFit<-eBayes(fit)

write.fit(ebFit, file=outfile , adjust="BH")
data<-read.table(outfile, sep="\t", header=T)

data<- topTable(ebFit, number=nrow(data))
write.csv(data,outfile)




