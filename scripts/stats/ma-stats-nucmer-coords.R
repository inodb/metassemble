#!/usr/bin/env Rscript
usage.string <-
"%prog [options] <assembled_contigs.fa> <mummer show-coords output file>

Displays some metagenomic related statistics after mapping assembled contigs
versus a set of references. The script is meant for determining the quality of
a metagenomic assembly of a known community. It uses the coords output of
show-coords from mummer, which should be used to map the contigs against the
reference.
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
if (length(args) != 2) {
    stop(paste(sep="\n","Two arguments required",usage.string))
}

coordstable <- read.table(args[[2]])
species <- unique(coordstable$V12)
reflengths <- unique(coordstable$V8)
names(reflengths) <- species

countgrep <- system(paste(sep="","grep -c '^>' ",args[[1]]),intern=TRUE)
count.coverage  <- function(coordstable, sp, qcovthreshold) {
    spcovsubset <- coordstable$V12 == sp & coordstable$V11 >= qcovthreshold
    differences <- coordstable$V1[spcovsubset][-1] -
        coordstable$V2[spcovsubset][-length(coordstable$V2[spcovsubset])]
    return(sum(coordstable$V5[spcovsubset]) +
           sum((differences - 1)[differences <= 0]))
}

# Contigs count table (contigs x reference)
# Number of times a contig is aligned to a reference
contigcounts <- table(coordstable$V13, coordstable$V12)
#                      dnn=list(contigs=coordstable$V13, species=coordstable$V12))
#maps.to.other.ref <- function(x, sp) {
#    return(length(x[!(names(x) %in% sp) & x > 0]))
#}
maps.to.multiple.ref <- function(x) {
    return(length(x[x > 0]))
}
map.ref.count <- apply(contigcounts, 1, maps.to.multiple.ref)

cat("Per genome\n")
cat(paste(sep="\t",
            "Reference",
            "1. Total number of bases",
            "2. Number of contigs mapping to this genome",
            "3. Nr of bases covered by (2) (Count overlapping bases once)",
            "4. Ratio (3 / 1)",
            "5. Nr of contigs with qcov 100%",
            "6. Nr of bases covered by (5) (Count overlapping bases once)",
            "7. Ratio (6 / 1)",
            "8. Nr of contigs mapping to multiple references",
            "9. Ratio (8 / 2)\n"))
for (sp in species) {
    qcov100 <- count.coverage(coordstable, sp, 100)
    qcovall <- count.coverage(coordstable, sp, 0)

    #browser()
    contigcountssp <- contigcounts[as.character(unique(coordstable$V13[coordstable$V12 == sp])),]
    # nr.of.other.ref.maps <- apply(contigcountssp, 1, maps.to.other.ref, sp)

    cat(sprintf("%s\t%i\t%i\t%i\t%.2f\t%i\t%i\t%.2f\t%i\t%.2f\n",
                  sp,
                  reflengths[sp],
                  nrow(contigcountssp),
                  qcovall,
                  qcovall / reflengths[sp] * 100.0,
                  length(unique(coordstable$V13[coordstable$V12 == sp &
                                coordstable$V11 == 100])),
                  qcov100,
                  qcov100 / reflengths[sp] * 100.0,
                  #length(nr.of.other.ref.maps[nr.of.other.ref.maps > 0]),
                  length(map.ref.count[names(map.ref.count) %in%
                    dimnames(contigcountssp)[[1]] & map.ref.count > 1]),
                  length(map.ref.count[names(map.ref.count) %in%
                    dimnames(contigcountssp)[[1]] & map.ref.count > 1]) /
                  nrow(contigcountssp) * 100.0 ))
}
cat("All genomes\n")
cat(paste(sep="\t",
            "Total number of contigs",
            "Total number of contigs mapping",
            "Total number of contigs mapping to one genome",
            "Total number of contigs mapping to multiple genomes\n"))
cat(sprintf("%i\t%i\t%i\t%i\n",
              as.integer(countgrep),
              nrow(contigcounts),
              length(map.ref.count[map.ref.count == 1]),
              length(map.ref.count[map.ref.count > 1])))
