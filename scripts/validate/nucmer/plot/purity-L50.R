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
if (length(args) < 2) {
    stop(paste(sep="\n","Number of arguments should be 2",usage.string))
}
output.pdf <- args[[1]]

assemblies <- read.table(args[[2]],header=TRUE)
assemblies$name <- sapply(as.character(assemblies$name), maxstring, 50)

print(levels(assemblies$asm_type))
print(assemblies$kmer_type[match(levels(assemblies$asm_type),assemblies$asm_type)])
print(as.numeric(assemblies$kmer_type[match(sort(levels(assemblies$asm_type)),assemblies$asm_type)]))

# Plot the assembly
pdf(output.pdf)
if (length(opts[['gradient-color-column']]) == 1) {
    g <- ggplot(data=assemblies,aes(x=as.numeric(l50),y=as.numeric(purity),size=as.integer(assemblies[,opts[['gradient-color-column']]]),color=asm_type,shape=asm_type))
    g <- g + geom_point() 
    g <- g + scale_size_continuous(name=opts[['gradient-color-column']]) + scale_colour_discrete("Assembly type") + scale_shape_manual("Assembly type",values=c(1,20)[as.numeric(assemblies$kmer_type[match(sort(levels(assemblies$asm_type)),assemblies$asm_type)])])
} else {
    g <- ggplot(data=assemblies,aes(x=as.numeric(l50),y=as.numeric(purity),color=asm_type))
    g <- g + geom_point() + scale_colour_discrete("Assembly")
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
