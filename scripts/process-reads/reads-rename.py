#! /usr/bin/env python
# Biopython version of
# https://github.com/ctb/khmer/blob/master/sandbox/multi-rename.py
"""Rename accessions. New names go from '>[prefix:]1[:original name]' to
'>[prefix:]n[:original name]' where n is the number of fasta sequences, [prefix:]
an optional prefix and [:original name] the optional original accession.

Usage:
    reads-rename.py [options] <fastafile>
Options:
    -c/--cutoff INT     Specify the cutoff [default 0].
    -h/--help           This help doc.
    -p/--prefix STRING  Prefix the accessions.
    -o/--original       Postfix with original name.
"""
import sys
import getopt
from Bio import SeqIO


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


class Error(Exception):
    def __init__(self, msg):
        self.msg = msg


def process(fastafile, prefix, cutoff, original):
    n = 0
    if prefix != '':
        prefix += ':'

    for record in SeqIO.parse(fastafile, "fasta"):
        if len(record.seq) >= cutoff:
            n += 1
            if original:
                print '>%s%s:%s\n%s' % (prefix, n, record.id, record.seq)
            else:
                print '>%s%s\n%s' % (prefix, n, record.seq)
 
    if n == 0:
        raise Error('No contigs >= %i\n' % cutoff)

def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hc:p:o", ["help", "cutoff",
                "prefix", "original"])
        except getopt.error, msg:
             raise Usage(msg)
        # process options
        cutoff = 0
        prefix = ''
        original = False
        for o, a in opts:
            if o in ("-h", "--help"):
                print __doc__
                return 0
            elif o in ("-c", "--cutoff"):
                cutoff = int(a)
            elif o in ("-p", "--prefix"):
                prefix = a
            elif o in ("-o", "--original"):
                original = True
        # process arguments
        if len(args) >= 1:
            for arg in args:
                process(fastafile=arg,
                        prefix=prefix,
                        cutoff=cutoff,
                        original=original)
        else:
            raise Usage("At least one fastafile required")
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, "for help use --help"
        return 2
    except Error, err:
        print >>sys.stderr, err.msg
        return 1

if __name__ == "__main__":
    sys.exit(main())
