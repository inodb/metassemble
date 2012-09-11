#!/usr/bin/env Rscript
usage.string <-
"%prog [options] <output.pdf> <ma-stats-nucmer-coords-purity.output>

Plots the output of ma-stats-nucmer-coords-purity in a plot.
"
suppressMessages({
    library(optparse)
    library(ggplot2)
})

maxstring <- function(string, cutoff) {
    if (nchar(string) > cutoff) {
        paste(sep="","...",substr(string,nchar(string)-cutoff+1,nchar(string)))
    } else {
        string
    }
}

# Parse arguments
parse_arguments <- function() {
  option_list <- list(make_option(c("-g", "--gradient-color-column"), type="character",
            default=NULL,
            help="Parses given column as a gradient to color the points by instead of discretely coloring the points by name."),
            make_option(c("-t","--title"), type="character",
            default="Purity-L50"))
  option_parser <- OptionParser(option_list=option_list,
                     usage=usage.string)
  opt <- parse_args(option_parser, positional_arguments=TRUE)
  return(opt)
}

opt <- parse_arguments()
args <- opt$args
opts <- opt$options
print(opts)
if (length(args) < 2) {
    stop(paste(sep="\n","Number of arguments should be 2",usage.string))
}
output.pdf <- args[[1]]

assemblies <- read.table(args[[2]],header=TRUE)
print(assemblies)
assemblies$name <- sapply(as.character(assemblies$name), maxstring, 50)

# Plot the assembly
pdf(output.pdf)
g <- ggplot(data=assemblies)
if (length(opts[['gradient-color-column']]) == 1) {
    g <- g + geom_point(aes(x=as.numeric(l50),y=as.numeric(purity),size=as.integer(assemblies[,opts[['gradient-color-column']]]),color=name,shape=kmer_type)) + 
        scale_size_continuous(name=opts[['gradient-color-column']]) + scale_colour_discrete("name") + scale_shape_manual(values=c(1,20))
} else {
    g <- g + geom_point(aes(x=as.numeric(l50),y=as.numeric(purity),color=name)) + 
        scale_colour_discrete("Assembly")
}
g <- g + 
    ylim(0,1) + 
    xlab("L50") +  
    ylab("Purity") +
    opts(legend.text=theme_text(size=8)) +
    opts(legend.direction="vertical",legend.position="right") +
    opts(title=opts[['title']])
print(g)
graphics.off()
