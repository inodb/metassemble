#!/usr/bin/env python
"""Newbler only handles contigs shorter than 2 kb. To do a combined assembly
(i.e. merge assemblies) one can split up contigs longer than 2 kb in shorter
overlapping ones and run Newbler on those. Something similar is done in [1].

First and last chunks of contigs > l are outputted twice, because Newbler
requires at least two reads to generate a consensus sequence. If those or not
outputted the beginning and end of the original contig are thus cut off when
re-assembling.

Usage:
cut-up-fasta.py <contigs.fasta> [contigs.fasta ...]
Options:
    -l INT   Cut contigs up in chunks of given parameter [1999]
    -o INT   Overlap between successive chunks. The last o bases of a chunk
             overlap with the first o bases of the next chunk. Has to be
             smaller than l. [1900]
Citations:
   1. Luo C, Tsementzi D, Kyrpides N, Read T, Konstantinidis KT (2012) Direct
   Comparisons of Illumina vs. Roche 454 Sequencing Technologies on the Same
   Microbial Community DNA Sample. PLoS ONE 7(2): e30087.
   doi:10.1371/journal.pone.0030087
"""
import sys
import getopt
from Bio import SeqIO

def cut_up_fasta(fastfiles, chunk_size, overlap):
    for ff in fastfiles:
        for record in SeqIO.parse(ff, "fasta"):
            if (len(record.seq) > chunk_size):
                i = 0
                # Output first part of contig twice
                print ">%s.%i\n%s" % (record.id, i,
                        record.seq[0:chunk_size])
                for split_seq in chunks(record.seq, chunk_size, overlap):
                    print ">%s.%i\n%s" % (record.id, i, split_seq)
                    i = i + 1
                # Output last chunk twice
                print ">%s.%i\n%s" % (record.id, i,
                        record.seq[-chunk_size + overlap:])
            else:
                print ">%s\n%s" % (record.id, record.seq)

    return 0

def chunks(l, n, o):
    """ Yield successive n-sized chunks from l with given overlap o between the
    chunks.
    """
    assert n > o

    for i in xrange(0, len(l) - o, n - o):
        yield l[i:i+n]

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hl:o:", ["help"])
        except getopt.error, msg:
             raise Usage(msg)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, "for help use --help"
        return 2
    # process options
    chunk_size = 1999
    overlap = 1900
    for o, a in opts:
        if o in ("-h", "--help"):
            print __doc__
            return 0
        if o in ("-l"):
            chunk_size = int(a)
        if o in ("-o"):
            overlap = int(a)
    if overlap >= chunk_size:
        print >>sys.stderr, "Overlap not smaller than chunk size"
        return 2
    # process arguments
    if (len(args) >= 1):
        return cut_up_fasta(args, chunk_size, overlap)
    else:
        print >>sys.stderr, "At least one argument required"
        print __doc__
        return 2

if __name__ == "__main__":
    sys.exit(main())
