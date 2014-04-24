import sys

from Bio import SeqIO
from Bio.SeqUtils import GC

# Argument one is a fasta file
for rec in SeqIO.parse(sys.argv[1], "fasta"):
    print "%s\t%1.2f\t%d" % (rec.id, GC(rec.seq) / 100.0, len(rec.seq))
