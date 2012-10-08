#!/usr/bin/env Rscript
usage.string <-
"%prog [options] <output.pdf> <ma-stats-nucmer-coords-purity.output>

Plots the output of ma-stats-nucmer-coords-purity in a plot.
"
suppressMessages({
    library(optparse)
    library(ggplot2)
})

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

# Plot the assembly
pdf(output.pdf)
if (length(opts[['gradient-color-column']]) == 1) {
    asm_type_hardcoded_levels <- 
        c('velvetnoscaf',
          'velvetscaf',
          'minimus2velvetnoscaf',
          'newblervelvetnoscaf',
          'bambus2velvetnoscaf',
          'bambus2minimus2velvetnoscaf',
          'bambus2newblervelvetnoscaf',
          'metavelvetnoscaf',
          'metavelvetscaf',
          'minimus2metavelvetnoscaf',
          'newblermetavelvetnoscaf',
          'bambus2metavelvetnoscaf',
          'bambus2minimus2metavelvetnoscaf',
          'bambus2newblermetavelvetnoscaf')
    assemblies$asm_type <- factor(assemblies$asm_type, levels=asm_type_hardcoded_levels)

    # Velvet filled circles, metavelvet unfilled circles
    shapes <- c(rep(19,7),rep(1,7))
    names(shapes) <- asm_type_hardcoded_levels
    # Color blind scale from:
    # http://wiki.stdout.org/rcookbook/Graphs/Colors%20%28ggplot2%29/
    # http://jfly.iam.u-tokyo.ac.jp/color/
    # noscaf black, scaf orange, minimus2noscaf light blue, newblernoscaf green,
    # bambus2noscaf dark blue, bambus2minimus red, bambus2newbler pink
    colors <- rep(c("#000000", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7"),2)
    names(colors) <- names(shapes)
    #colorfill <- c("#E69F00","#CC79A7","#000000","#0072B2",NA,NA,NA,NA)
    #names(colorfill) <- names(shapes)

    g <- ggplot(data=assemblies,aes(x=as.numeric(l50),y=as.numeric(purity),size=assemblies[,opts[['gradient-color-column']]],color=asm_type,shape=asm_type))
    # scale_fill_manual("Assembly",values=colorfill[levels(assemblies$asm_type)]) +
    g <- g + geom_point() +
        scale_size_continuous(name=opts[['gradient-color-column']],range=c(1,6)) +
        scale_colour_manual("Assembly",values=colors[levels(assemblies$asm_type)]) +
        scale_shape_manual("Assembly",values=shapes[levels(assemblies$asm_type)])
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
