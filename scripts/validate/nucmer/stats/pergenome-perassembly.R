#!/usr/bin/env Rscript
usage.string <-
"%prog [options] <assembled_contigs.fa>
                 <references.stats>
                 <mummer show-coords output file>
                 <output dir>

Displays some metagenomic related statistics after mapping assembled contigs
versus a set of references. The script is meant for determining the quality of
a metagenomic assembly of a known community. It uses the coords output of
show-coords from mummer, which should be used to map the contigs against the
reference.
"
suppressMessages({
    library(optparse)
})

computeNonOverlappingBasesQueryFwd <- function(df) {
    # Computes non overlapping bases per contig per genome in the forward
    # direction to determine purity of a contig. Use in combination with ddply.
    # Reverse direction can be computed by simply swappig columns V4 and V3.
    # There are numerous problems with this approach, so using the longest
    # alignment seems like a more sensible choice.
    # TODO: It only looks at overlap in the contig, not the reference genome,
    # so a tandem repeat in the contig, with a single occurence in the
    # reference genome gives 100% purity over that area.
    # TODO: Which alignment should be used with overlapping alignment? Now the
    # the alignments are ordered on lower bound and the lower bounds of all
    # following rows are increased by the amount of overlap i.e. the overlap is
    # removed from each subsequent entry.
    # TODO: When parts of the contig map to totally different locations it is
    # still considered pure while that does not seem very logical.
    # TODO: This is extremely slooooooooow.
    if (nrow(df.fwd) == 0) {
        return(0) 
    }
    df.fwd <- subset(df, V4 > V3)
    if (nrow(df.fwd) == 0) {
        return(0) 
    }
    df.fwd <- df.fwd[order(df.fwd$V3),]
    repeat {
        if (nrow(df.fwd) == 1) {
            return(df.fwd$V11) 
        }
        contained.alignments.index <- which(df.fwd$V4[-nrow(df.fwd)] -
            df.fwd$V4[-1] >= 0) + 1
        if (length(contained.alignments.index) > 0) {
            df.fwd <- df.fwd[-contained.alignments.index,] 
        } else {
            break 
        } 
    }

    overlap <- c(0, sapply(df.fwd$V4[-nrow(df.fwd)] - df.fwd$V3[-1] + 1, 
        max, 0))
    sum((df.fwd$V6 - overlap) * df.fwd$V7 / df.fwd$V9)  
}


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
if (length(args) != 4) {
    stop(paste(sep="\n","Four arguments required",usage.string))
}

coordstable <- read.table(args[[3]])
species <- unique(coordstable$V12)
#reflengths <- unique(coordstable$V8)
#names(reflengths) <- species

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

ref.stats <- read.table(args[[2]],header=TRUE)

pergenome.stats <- data.frame(matrix(nrow=0,ncol=11))
sink(paste(sep="",args[[4]],'/gnm-stats.tab'))
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
            "9. Ratio (8 / 2)",
            "10. GC-content of reference",
            "11. Original read coverage\n"))
for (sp in species) {
    qcov100 <- count.coverage(coordstable, sp, 100)
    qcovall <- count.coverage(coordstable, sp, 0)

    contigcountssp <- contigcounts[as.character(unique(coordstable$V13[coordstable$V12 == sp])),]

    cat(sprintf("%s\t%i\t%i\t%i\t%.2f\t%i\t%i\t%.2f\t%i\t%.2f\t%.2f\t%.2f\n",
                  sp,
                  ref.stats[sp,]$length,
                  nrow(contigcountssp),
                  qcovall,
                  qcovall / ref.stats[sp,]$length,
                  length(unique(coordstable$V13[coordstable$V12 == sp &
                                coordstable$V11 == 100])),
                  qcov100,
                  qcov100 / ref.stats[sp,]$length,
                  length(map.ref.count[names(map.ref.count) %in%
                    dimnames(contigcountssp)[[1]] & map.ref.count > 1]),
                  length(map.ref.count[names(map.ref.count) %in%
                    dimnames(contigcountssp)[[1]] & map.ref.count > 1]) /
                  nrow(contigcountssp),
                  ref.stats[sp,]$GC_content,
                  ref.stats[sp,]$cov))
}
sink()
# TODO: Per genome stats are read in here with the idea that the entire
# previous part should be computed in some other way that is less horrendous.
pergenome.stats <- read.table(paste(sep="",args[[4]],'/gnm-stats.tab'),skip=1)

sink(paste(sep='',args[[4]],'/asm-stats.tab'))
# cat(file=args[[4]],"All genomes\n")
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
sink()

# GC content versus contig coverage plot
sink("/dev/null")
pdf(paste(sep='',args[[4]],'/gccontent-contigcov.pdf'))
plot(pergenome.stats$V10,pergenome.stats$V7,xlim=c(0,1),ylim=c(0,1),xlab="GC content",ylab="Contig coverage")
dev.off()

# Purity-Length Plot
max.qcov.per.contig <- aggregate(coordstable$V11, list(coordstable$V13), max)
names(max.qcov.per.contig) <- c('name','max')
max.qcov.per.contig$contigLength <- coordstable$V9[match(max.qcov.per.contig$name, coordstable$V13)]
pdf(paste(sep='',args[[4]],'/purity-length.pdf'))
length.breaks <- hist(max.qcov.per.contig$contigLength, plot=FALSE)
length.bounds.matrix <- cbind(
    lower=length.breaks$breaks[-length(length.breaks$breaks)],
    upper=length.breaks$breaks[-1])
mean.max.contig.cov.by.length <- apply(length.bounds.matrix, 1, function(x) {
    colMeans(subset(max.qcov.per.contig,
        contigLength > x['lower'] & contigLength < x['upper'],
        select=max))
})
write.table(data.frame(lower=length.bounds.matrix[,'lower'],mid=length.breaks$mid,upper=length.bounds.matrix[,'upper'],counts=length.breaks$counts,mean_max_contig_cov=mean.max.contig.cov.by.length),quote=FALSE,file=paste(sep='',args[[4]],'/purity-length-hist.tab'))
plot(length.breaks$mid,mean.max.contig.cov.by.length,lwd=log(length.breaks$counts),ylab="Purity",xlab="Contig length",ylim=c(0,100))
dev.off()
