#!/usr/bin/env Rscript
usage.string <- 
"%prog [options] <contig-length-cut-off> <contigs.fa> [contigs2.fa ...]
or
%prog [options] <contig-length-cut-off> <label>:<contigs.fa> [label2:contigs2.fa ...]
Plots histogram of contig lengths. If multiple fasta files are given the
combined dataset is used to compute the number of bins using Sturges' formula.
Each histogram is then plotted on top of each other using the calculated number
of bins for all datapoints combined.
"

# Parse arguments
parse_arguments <- function() {
  suppressMessages({
    library(optparse)
    library(ShortRead)
    library(ggplot2)
  })
  option_list <- list()
  option_parser <- OptionParser(option_list=option_list,
                     usage=usage.string)
  opt <- parse_args(option_parser, positional_arguments=TRUE)
  return(opt)
}

maxstring <- function(string, cutoff) {
    if (nchar(string) > cutoff) {
        paste(sep="","...",substr(string,nchar(string)-cutoff+1,nchar(string)))
    } else {
        string
    }
}

opt <- parse_arguments()
args <- opt$args
if (length(args) == 0) {
    stop(paste(sep="\n","Zero arguments",usage.string))
}

# Extract label:file combinations from arguments
files <- c()
labels <- c()
cutoff <- as.integer(args[[1]])
# From is.integer example
is.wholenumber <-
        function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
if (!is.wholenumber(cutoff)) {
    stop(paste("Given cut off", cutoff,
               "is not a whole number.\n",
               "See -h for help"))
}
for (arg in args[2:length(args)]) {
    splitv <-  strsplit(arg,":")[[1]]
    if (length(splitv) == 2) {
        files <- c(files, splitv[[2]])
        labels <- c(labels,splitv[[1]])
    } else if (length(splitv) == 1) {
        files <- c(files, splitv[[1]])
        labels <- c(labels, maxstring(splitv[[1]],50))
    } else {
        stop(paste("Argument", splitlist,
                   "is not of type label:file.\n",
                   "See -h for help"))
    }
}

alllengths <- c()
nrofsequences <- c()
for (file in files) {
    tmp <- fasta.info(file)
    alllengths <- c(alllengths, tmp[tmp>cutoff])
    nrofsequences <- c(nrofsequences, length(tmp[tmp>cutoff]))  
}

Data <- data.frame(assembly=rep(labels,nrofsequences),contiglength=alllengths)
histall <- hist(alllengths,plot=FALSE)
# Also include upper boundary in the breaks
breaksall <- c(histall$breaks,                           
               histall$breaks[length(histall$breaks)] +  
               histall$breaks[2]-histall$breaks[1])
hists <- tapply(Data$contiglength, Data$assembly,
                function(i) {
                    tmp <- hist(i,
                                breaks=breaksall,
                                plot=FALSE)
                    data.frame(br=tmp$breaks,co=c(tmp$counts,0))
                })
print(hists)
ll <- sapply(hists,nrow)
hists <- do.call(rbind, hists)
browser()
hists$fac <- rep(as.character(labels),ll)
png("test.png")
print(qplot(br,co,data=hists,geom="step",colour=as.character(fac),
            xlab="Contig length",ylab="Total bases in contigs with length x")
      +
      scale_colour_discrete(name="Contigs") +
      opts(legend.direction="vertical"))#legend.position="bottom"
graphics.off()
browseURL("test.png")
