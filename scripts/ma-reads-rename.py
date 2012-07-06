#! /usr/bin/env python
# Biopython version of
# https://github.com/ctb/khmer/blob/master/sandbox/multi-rename.py
import sys
from Bio import SeqIO

CUTOFF=200

n = 0
for filename in sys.argv[1:]:
   for record in SeqIO.parse(filename, "fasta"):
       if len(record.seq) >= CUTOFF:
           n += 1
           print '>%s %s\n%s' % (n, record.id, record.seq)

if n == 0:
    sys.stderr.write('ERROR: No contigs >= %i\n' % CUTOFF)
    sys.exit(1)
