import argparse
import sqlite3

COORDS = "/bubo/glob/g16/inod/metagenomics/result/chris-mock/Project_ID793_dAmore/Sample_500pg_unbalanced/ma2-out/assemblies/velvet/noscaf/noscaf_31/val/nucmer.coords"
REFSTATS = "/bubo/glob/g16/inod/metagenomics/result/chris-mock/Project_ID793_dAmore/Sample_500pg_unbalanced/ma2-out/reference-stats/ref.stats"


def calc_genome_contig_cov_in_bases(cursor):
    """Genome contig coverage is a metric that indicates how well the genome is
    covered by contigs. For each contig only the purest alignment is
    considered. Purity is defined as COVQ * IDY / 10,000. If there are muliple
    alignments for a contig with maximum purity, both are used. Covered bases
    are only counted once and computed by multiplying the length of the
    alignment in the reference (S1 - E1 + 1) by the alignment identity (IDY).
    If contigs overlap the first contig based on its location in the reference
    genome is used. If the second contig extends further than the first, the
    second's bases are added as well using the IDY of the second. One might
    prefer to use the purest alignment in case of overlap but I'll implement
    that only if persuaded with dinner.

    Returns a dictionary of the genome contig coverage per genome as the number
    of non-overlapping bases that are covered by contigs aligning with maximum
    purity.
    """
    refcov = {}  # nr of bases covered by contigs aligned with max purity
    prev_e1 = 0
    for row in cursor.execute("""Select * From ( Select *, max(COVQ * IDY /
                              10000) As purity FROM Coords Group by QRYID)
                              Where COVQ * IDY / 10000 == purity Order by
                              REFID, S1 asc"""):
        if row["REFID"] in refcov:
            if prev_e1 >= row["E1"]:
                continue
            elif prev_e1 >= row["S1"]:
                refcov[row["REFID"]] += int((row["E1"] - prev_e1) * row["IDY"])
            else:
                refcov[row["REFID"]] += int((row["E1"] - row["S1"] + 1) *
                                            row["IDY"])
        else:
            refcov[row["REFID"]] = int((row["E1"] - row["S1"]) * row["IDY"]) + 1
        prev_e1 = row["E1"]

    return(refcov)


def readtable(tablefile, sep=None):
    """Reads given table separated 'sep'.

    Returns a two-dimensional dictionary. Outer dictionary uses first column in
    'tablefile' as key, inner dictionary uses column names found in the first
    line as column names. The number of column names should be equal to number
    of columns - 1."""
    table2d = {}

    # Get column names
    tfh = open(tablefile, "r")
    line = tfh.readline()
    cols = line.split()

    # Insert rows
    for line in tfh:
        splits = line.split()
        table2d[splits[0]] = {}
        for i in range(1, len(splits)):
            table2d[splits[0]][cols[i - 1]] = splits[i]

    return(table2d)


def main(coordsfile, refstatsfile):
    # Create db
    #dbc = sqlite3.connect(':memory:')
    dbc = sqlite3.connect(coordsfile + ".sqlite")
    try:
        dbc.execute("""Create table Coords (ID INTEGER PRIMARY KEY, S1 INTEGER, E1
                    INTEGER, S2 INTEGER, E2 INTEGER, LEN1 INTEGER, LEN2 INTEGER,
                    IDY REAL, LENR INTEGER, LENQ INTEGER, COVR REAL, COVQ REAL,
                    REFID TEXT, QRYID TEXT)""")
        dbc.row_factory = sqlite3.Row

        # Parse file and add to db
        columns = ('S1', 'E1', 'S2', 'E2', 'LEN1', 'LEN2', 'IDY', 'LENR',
                   'LENQ', 'COVR', 'COVQ', 'REFID', 'QRYID')
        cfh = open(coordsfile, "r")
        for line in cfh:
            parsed_line = dict(zip(columns, line.split()))
            dbc.execute("""Insert into Coords values(NULL, :S1, :E1, :S2, :E2,
                        :LEN1, :LEN2, :IDY, :LENR, :LENQ, :COVR, :COVQ, :REFID,
                        :QRYID)""", parsed_line)
    except:
        dbc.row_factory = sqlite3.Row
        pass

    gconcov = calc_genome_contig_cov_in_bases(dbc.cursor())
    dbc.commit()

    # Get all reference lengths
    reflens = readtable(refstatsfile)

    # Add 0 covered bases for genomes with no max purity aligned contigs
    for genome in reflens.keys():
        if genome not in gconcov.keys():
            gconcov[genome] = 0

    # Print genome contig coverage table
    print("Genome_contig_cov_bases Genome_length Genome_contig_cov_ratio")
    for genome in reflens:
        print("%s %i %i %f" % (genome, gconcov[genome],
                               int(reflens[genome]["length"]),
                               gconcov[genome] /
                               int(reflens[genome]["length"])))
    sum_cov_bases = sum(gconcov.values())
    sum_bases = sum(int(ref["length"]) for ref in reflens.values())
    print("%s %i %i %f" % ("Total", sum_cov_bases, sum_bases, sum_cov_bases /
                           sum_bases))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("coordsfile", help="Output from nucmer \
        show-coords -h")
    parser.add_argument("refstatsfile", help="Reference stats file")
    #args = parser.parse_args()
    #main(args.coordsfile, args.refstatsfile)
    main(COORDS, REFSTATS)
