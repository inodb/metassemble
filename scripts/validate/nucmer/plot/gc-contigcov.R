#!/usr/bin/env Rscript
usage.string <-
"%prog [options] <per-genome-stats.out> 
"
suppressMessages({
    library(optparse)
})

# Parse arguments
parse_arguments <- function() {
  option_list <- list()
  option_parser <- OptionParser(option_list=option_list,
                     usage=usage.string)
  opt <- parse_args(option_parser, positional_arguments=TRUE)
  return(opt)
}

opt <- parse_arguments()
args <- opt$args
if (length(args) != 1) {
    stop(paste(sep="\n","One argument required",usage.string))
}
