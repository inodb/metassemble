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
    # Velvet circles, metavelvet diamonds
    shapes <- c(20,20,20,20,18,18,18,18)
    names(shapes) <- c('velvetscaf',
        'newblervelvetnoscaf','velvetnoscaf','minimus2velvetnoscaf',
        'metavelvetscaf','metavelvetnoscaf','minimus2metavelvetnoscaf','newblermetavelvetnoscaf')
    # Color blind scale from:
    # http://wiki.stdout.org/rcookbook/Graphs/Colors%20%28ggplot2%29/
    # http://jfly.iam.u-tokyo.ac.jp/color/
    # noscaf black, scaf yellowish, minimus2noscaf blue, newblernoscaf pink 
    colors <- c("#E69F00","#CC79A7","#000000","#0072B2","#E69F00","#000000","#0072B2","#CC79A7")
    names(colors) <- names(shapes)

    g <- ggplot(data=assemblies,aes(x=as.numeric(l50),y=as.numeric(purity),size=assemblies[,opts[['gradient-color-column']]],color=asm_type,shape=asm_type,fill=asm_type))
    g <- g + geom_point() +
        scale_size_continuous(name=opts[['gradient-color-column']],range=c(1,6)) +
        scale_colour_manual("Assembly",values=colors[levels(assemblies$asm_type)]) +
        scale_fill_manual("Assembly",values=colors[levels(assemblies$asm_type)]) +
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
