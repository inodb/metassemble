import argparse
import sqlite3
import re
import os
import errno
import pysam
import matplotlib.pyplot as plt
import numpy
from collections import Counter
from collections import defaultdict  # is faster than Counter

#CONTIGS = "/bubo/glob/g16/inod/metagenomics/result/chris-mock/Project_ID793_dAmore/Sample_500pg_unbalanced/ma2-out/assemblies/velvet/noscaf/noscaf_31/contigs.fa"
#COORDS = "/bubo/glob/g16/inod/metagenomics/result/chris-mock/Project_ID793_dAmore/Sample_500pg_unbalanced/ma2-out/assemblies/velvet/noscaf/noscaf_31/val/nucmer.coords"
#REFSTATS = "/bubo/glob/g16/inod/metagenomics/result/chris-mock/Project_ID793_dAmore/Sample_500pg_unbalanced/ma2-out/reference-stats/ref.stats"
#OUTDIR = "/bubo/glob/g16/inod/metagenomics/result/chris-mock/Project_ID793_dAmore/Sample_500pg_unbalanced/ma2-out/assemblies/velvet/noscaf/noscaf_31/val/valnew"
#BAMREF = "/bubo/glob/g16/inod/metagenomics/result/chris-mock/Project_ID793_dAmore/Sample_500pg_unbalanced/ma2-out/reference-stats/ref_500pg_unbalanced.qtrim-smds.bam"
#BAMASM = "/bubo/glob/g16/inod/metagenomics/result/chris-mock/Project_ID793_dAmore/Sample_500pg_unbalanced/ma2-out/assemblies/velvet/noscaf/noscaf_31/bambus2/contigs_500pg_unbalanced-smds.bam"


def make_dir(directory):
    """Make directory unless existing. Ignore error in the latter case."""
    try:
        os.makedirs(directory)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def calc_genome_contig_cov_in_bases(cursor, cut_off=100):
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
    for row in cursor.execute("""SELECT * FROM (SELECT *, max(COVQ * IDY /
                              10000) AS purity FROM Coords GROUP BY QRYID)
                              WHERE COVQ * IDY / 10000 == purity AND LENQ >= :cut_off ORDER BY
                              REFID, S1 ASC""", dict(cut_off=cut_off)):
        if row["REFID"] in refcov:
            if prev_e1 >= row["E1"]:
                continue
            elif prev_e1 >= row["S1"]:
                refcov[row["REFID"]] += int((row["E1"] -
                                            prev_e1) * (row["IDY"] / 100.0))
            else:
                refcov[row["REFID"]] += int((row["E1"] - row["S1"] + 1) *
                                            (row["IDY"] / 100.0))
        else:
            refcov[row["REFID"]] = int((row["E1"] - row["S1"])
                                       * (row["IDY"] / 100.0)) + 1
        prev_e1 = row["E1"]

    return(refcov)


def calc_max_purity_per_contig(cursor, cut_off=100):
    q_max_purity = {}
    for row in cursor.execute("""SELECT * FROM (SELECT *, max(COVQ * IDY /
                              10000) AS purity FROM Coords GROUP BY QRYID)
                              WHERE COVQ * IDY / 10000 == purity
                              AND LENQ >= :cut_off""", dict(cut_off=cut_off)):
        q_max_purity[row["QRYID"]] = dict(max_purity=row["purity"],
                                          length=row["LENQ"])

    return(q_max_purity)


def calc_alignedbases_per_contig(cursor, cut_off=100):
    q_aln_bases = {}  # nr of bases aligned per contig
    prev_end = 0
    for row in cursor.execute("""SELECT *, min(S2, E2) AS start, max(S2, E2) as end FROM COORDS WHERE LENQ >= :cut_off ORDER BY
                              QRYID, start ASC""", dict(cut_off=cut_off)):
        if row["QRYID"] in q_aln_bases:
            if prev_end >= row["end"]:
                continue
            elif prev_end >= row["start"]:
                q_aln_bases[row["QRYID"]] += row["end"] - prev_end
            else:
                q_aln_bases[row["QRYID"]] += row["end"] - row["start"] + 1
        else:
            q_aln_bases[row["QRYID"]] = row["end"] - row["start"]
        prev_end = row["end"]

    return(q_aln_bases)


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
    cols = line.strip().split(sep)

    # Insert rows
    for line in tfh:
        splits = line.strip().split(sep)
        table2d[splits[0]] = {}
        for i in range(1, len(splits)):
            table2d[splits[0]][cols[i - 1]] = splits[i]

    return(table2d)


def asm_stats_fasta(fafile, cut_off=100):
    """Return L50 of an assembly and the total number of bases. 50% of all
    bases in the assembly are located in contigs equal or larger than l50."""
    name = None
    all_contig_lengths = {}

    # Determine lengths of all contigs
    for line in open(fafile):
        if line[0] == '>':
            name = line[1:-1]
            all_contig_lengths[name] = 0
        else:
            all_contig_lengths[name] += len("".join(re.findall(
                "[ACGTacgt]+", line)))

    # Determine l50, 50% of the bases in the fasta file are located in contigs => l50
    sorted_contigs = sorted(filter(lambda x: x >= cut_off,
                                   all_contig_lengths.values()),
                            reverse=True)
    cumsum = 0
    totbases = sum(sorted_contigs)
    max_length = sorted_contigs[0]
    for i in range(len(sorted_contigs)):
        cumsum += sorted_contigs[i]
        if cumsum >= totbases / 2:
            l50 = sorted_contigs[i]
            n50 = i
            break

    return(l50, n50, totbases, max_length, len(sorted_contigs), all_contig_lengths)


def calc_read_level_purity(bamref, bamasm, refphylfile, cut_off=100):
    reffile = pysam.Samfile(bamref, "rb")
    # Determine which reads map uniquely to what reference
    read_ref_map = {}
    for record in reffile:
        if record.mapq > 0:
            if record.qname in read_ref_map:
                if record.tid not in read_ref_map[record.qname]:
                    read_ref_map[record.qname].append(record.tid)
            else:
                read_ref_map[record.qname] = [record.tid]

    # Determine for every contig (or scaffold) which reads map to what reference
    asmfile = pysam.Samfile(bamasm, "rb")
    contig_ref_map = {}
    for record in asmfile:
        if record.mapq > 0 and record.qlen >= cut_off:
            # The read should be mapped to the reference genomes
            if record.qname in read_ref_map:
                # Check unambigiuously mapping reads and ambiguously mapping reads
                if not record.tid in contig_ref_map:
                    contig_ref_map[record.tid] = {}
                    contig_ref_map[record.tid]["unamb"] = Counter()
                    contig_ref_map[record.tid]["amb"] = Counter()

                if len(read_ref_map[record.qname]) == 1:
                    ref_tid = read_ref_map[record.qname][0]
                    contig_ref_map[record.tid]["unamb"][reffile.getrname(ref_tid)] += 1
                    contig_ref_map[record.tid]["amb"][reffile.getrname(ref_tid)] += 1
                else:
                    for ref_tid in read_ref_map[record.qname]:
                        contig_ref_map[record.tid]["amb"][reffile.getrname(ref_tid)] += 1

    # Determine read level purity for every contig i.e. (number of reads
    # mapping to most dominant strain / number of reads mapping to contig).
    # Both for unambiguous and ambiguous reads. Also determine nr reads mapping
    # at each taxonomic level to give an indication at what level chimericity
    # occurs
    contig_read_level_purity = {}
    rp = numpy.genfromtxt(refphylfile, names=True, dtype=None,
                          missing_values="-", delimiter="\t")
    tax_lvls = [rp.dtype.names[i] for i in range(4, 14)]
    for k, v in contig_ref_map.iteritems():
        contig_read_level_purity[asmfile.getrname(k)] = {}

        # Calculate for both "unamb" and "amb"
        for kj, vj in v.iteritems():
            # Get most dominant strain and number of reads
            mc = vj.most_common(1)
            if len(mc) != 0:
                dom_nr_reads = mc[0][1]
                tot_nr_reads = sum(vj.itervalues())
                read_level_purity = float(dom_nr_reads) / tot_nr_reads
                dominant_strain = mc[0][0]

                contig_read_level_purity[asmfile.getrname(k)][kj] = \
                    dict(read_level_purity=read_level_purity,
                         dominant_strain=dominant_strain,
                         tot_nr_reads=tot_nr_reads,
                         dom_nr_reads=dom_nr_reads)

                # Determine nr of reads per taxon and the taxonomic level of
                # the lowest common ancestor (LCA)
                if dom_nr_reads != tot_nr_reads:
                    nr_reads_per_taxon, lca = calc_nr_reads_per_taxon(vj, rp,
                                                                      tax_lvls)
                    contig_read_level_purity[asmfile.getrname(k)][kj]["lca"] =\
                            lca
                else:
                    nr_reads_per_taxon = [(tl, tot_nr_reads) for tl in
                                          tax_lvls]
                    contig_read_level_purity[asmfile.getrname(k)][kj]["lca"] =\
                            tax_lvls[-1]
                contig_read_level_purity[asmfile.getrname(k)][kj].update(nr_reads_per_taxon)

    return(contig_read_level_purity)


def calc_nr_reads_per_taxon(refcounts, rp, tax_lvls):
    """Determines the number of reads per taxon if the reads map to multiple
    strains (otherwise the number would be equal for all taxa). Returns a tuple
    of two elements, first item is  a list of (taxon, number of reads) in
    decreasing taxonomic rank, the second the taxonomic rank of the lowest
    common ancestor (LCA) of the reads.

    Keyword arguments:
    refcounts -- dictionary with 'reference : count' pairs
    rp        -- numpy array of references phylogeny
    tax_lvls  -- list of taxonomic levels, must be in rp.dtype.names
    """
    mc = refcounts.most_common(1)
    assert(len(mc) == 1)

    dom_str = rp["fasta_name"] == mc[0][0]
    assert(len(rp[dom_str]) == 1)  # Strain should appear only once in array

    # Count number of reads at each taxon. Start from lowest taxonomic level to
    # highest and stop at sum
    nr_reads_per_taxon = []
    sum_reads = sum(refcounts.itervalues())
    taxon_sum_reads = -1
    for tl in reversed(tax_lvls):
        if taxon_sum_reads < sum_reads:
            dom_str_taxon = rp[dom_str][tl]
            taxon_sum_reads = sum(v for k, v in refcounts.iteritems() if
                                  rp[rp["fasta_name"] == k][tl] ==
                                  dom_str_taxon)
            lca = tl
        nr_reads_per_taxon.append((tl, taxon_sum_reads))

    # Reverse back so taxa or from highest to lowest rank
    nr_reads_per_taxon.reverse()

    return(nr_reads_per_taxon, lca)


def print_dict2tsv(d, filepath):
    with open(filepath, "w") as fh:
        fh.write("\t".join(d.keys()) + "\n")
        fh.write("\t".join(str(v) for v in d.itervalues()) + "\n")


def get_default_nested_dict(d, keys, default="-", map_values=False):
    di = d

    for k in keys:
        try:
            di = di[k]
        except KeyError:
            return default

    if map_values:
        return map_values.get(di, di)
    else:
        return di


def main(coordsfile, refstatsfile, refphylfile, contigs, bamref, bamasm,
         outdir, name="-", asm_type="-", kmer_type="-", kmer_size="-",
         kmin="-", kmax="-", cut_off=100):
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

    dbc.commit()
    cur = dbc.cursor()
    gconcov = calc_genome_contig_cov_in_bases(cur)
    q_aln_bases = calc_alignedbases_per_contig(cur)
    contig_max_purity = calc_max_purity_per_contig(cur)
    sum_purest_bases = sum(
        c["max_purity"] * c["length"] for c in contig_max_purity.itervalues())
    sum_aln_bases = sum(q_aln_bases.itervalues())

    l50, n50, q_sum_bases, q_max, trim_n, all_contig_lengths = asm_stats_fasta(
        contigs)
    # Alignment purity, how pure are all the alignments
    aln_pur = float(sum_purest_bases) / sum_aln_bases
    # Global purity
    glob_pur = float(sum_purest_bases) / q_sum_bases
    aln_ratio = float(sum_aln_bases) / q_sum_bases

    # Get all reference lengths
    reflens = readtable(refstatsfile, sep="\t")
    # Add 0 covered bases for genomes with no max purity aligned contigs
    for genome in reflens.keys():
        if genome not in gconcov.keys():
            gconcov[genome] = 0
    ref_sum_bases = sum(int(ref["length"]) for ref in reflens.values())

    make_dir(outdir)

    # Print assembly stats
    print_dict2tsv(d=dict(name=name, global_purity=glob_pur,
                          aln_purity=aln_pur, aln_ratio=aln_ratio, l50=l50,
                          n50=n50, trim_n=trim_n, max_contig_length=q_max,
                          cut_off=100, trim_n_mapping=len(q_aln_bases),
                          sum_ref_lengths=ref_sum_bases,
                          sum_purest_bases=sum_purest_bases,
                          metagenome_cov=sum_purest_bases / ref_sum_bases,
                          sum_bases=q_sum_bases, asm_type=asm_type,
                          kmer_type=kmer_type, kmer_size=kmer_size, kmin=kmin,
                          kmax=kmax), filepath=outdir + "/asm-stats.tsv")

    # Print genome contig coverage
    with open(outdir + "/genome-contig-coverage.tsv", "w") as fh:
        fh.write("genome\tgenome_contig_cov_bases\tgenome_length\tgenome_contig_cov_ratio\tGC_content\tread_cov_ratio\n")
        for genome in reflens:
            fh.write("%s\t%i\t%i\t%f\t%f\t%f\n" % (genome, gconcov[genome],
                     int(reflens[genome]["length"]), gconcov[genome] /
                     float(reflens[genome]["length"]),
                     float(reflens[genome]["GC_content"]),
                     float(reflens[genome]["ratio_covered"])))

    # Print contig purity
    contig_read_level_purity = calc_read_level_purity(bamref, bamasm,
                                                      refphylfile, cut_off)
    with open(outdir + "/contig-purity.tsv", "w") as fh:
        fh.write("contig\tunamb_read_level_purity\tunamb_tot_nr_reads\t"
                 "unamb_dom_nr_reads\t"
                 "unamb_dominant_strain\tamb_read_level_purity\t"
                 "amb_tot_nr_reads\tamb_dom_nr_reads\t"
                 "amb_dominant_strain\tmax_aln_purity\t"
                 "contig_length\t"
                 "unamb_nr_reads_life\tamb_nr_reads_superkingdom\t"
                 "unamb_nr_reads_superphylum\tamb_nr_reads_phylum\t"
                 "unamb_nr_reads_class\tamb_nr_reads_order\t"
                 "unamb_nr_reads_family\tamb_nr_reads_genus\t"
                 "unamb_nr_reads_species\tunamb_nr_reads_strain\t"
                 "unamb_lca\t"
                 "amb_nr_reads_life\tamb_nr_reads_superkingdom\t"
                 "amb_nr_reads_superphylum\tamb_nr_reads_phylum\t"
                 "amb_nr_reads_class\tamb_nr_reads_order\t"
                 "amb_nr_reads_family\tamb_nr_reads_genus\t"
                 "amb_nr_reads_species\tamb_nr_reads_strain\t"
                 "amb_lca"
                 "\n")
        for c in dict((k, v) for k, v in all_contig_lengths.iteritems() if v >= cut_off):
            fh.write(
                "\t".join(str(x) for x in
                [c,
                 get_default_nested_dict(contig_read_level_purity, [c, "unamb",
                                                                    "read_level_purity"]),
                 get_default_nested_dict(contig_read_level_purity, [c, "unamb",
                                                                    "tot_nr_reads"]),
                 get_default_nested_dict(contig_read_level_purity, [c, "unamb",
                                                                    "dom_nr_reads"]),
                 get_default_nested_dict(contig_read_level_purity, [c, "unamb",
                                                                    "dominant_strain"]),
                 get_default_nested_dict(contig_read_level_purity, [c, "amb",
                                                                    "read_level_purity"]),
                 get_default_nested_dict(contig_read_level_purity, [c, "amb",
                                                                    "tot_nr_reads"]),
                 get_default_nested_dict(contig_read_level_purity, [c, "amb",
                                                                    "dom_nr_reads"]),
                 get_default_nested_dict(contig_read_level_purity, [c, "amb",
                                                                    "dominant_strain"]),
                 get_default_nested_dict(contig_max_purity, [c, "max_purity"]),
                 all_contig_lengths[c],
                 get_default_nested_dict(contig_read_level_purity, [c, "unamb",
                                                                     "no_rank"]),
                 get_default_nested_dict(contig_read_level_purity, [c, "unamb",
                                                                     "superkingdom"]),
                 get_default_nested_dict(contig_read_level_purity, [c, "unamb",
                                                                     "superphylum"]),
                 get_default_nested_dict(contig_read_level_purity, [c, "unamb",
                                                                     "phylum"]),
                 get_default_nested_dict(contig_read_level_purity, [c, "unamb",
                                                                     "class"]),
                 get_default_nested_dict(contig_read_level_purity, [c, "unamb",
                                                                     "order"]),
                 get_default_nested_dict(contig_read_level_purity, [c, "unamb",
                                                                     "family"]),
                 get_default_nested_dict(contig_read_level_purity, [c, "unamb",
                                                                     "genus"]),
                 get_default_nested_dict(contig_read_level_purity, [c, "unamb",
                                                                     "species"]),
                 get_default_nested_dict(contig_read_level_purity, [c, "unamb",
                                                                     "topname"]),
                 get_default_nested_dict(contig_read_level_purity, [c, "unamb",
                                                                    "lca"],
                                         map_values=dict(topname="strain",
                                                         no_rank="life")),
                 get_default_nested_dict(contig_read_level_purity, [c, "amb",
                                                                     "no_rank"]),
                 get_default_nested_dict(contig_read_level_purity, [c, "amb",
                                                                     "superkingdom"]),
                 get_default_nested_dict(contig_read_level_purity, [c, "amb",
                                                                     "superphylum"]),
                 get_default_nested_dict(contig_read_level_purity, [c, "amb",
                                                                     "phylum"]),
                 get_default_nested_dict(contig_read_level_purity, [c, "amb",
                                                                     "class"]),
                 get_default_nested_dict(contig_read_level_purity, [c, "amb",
                                                                     "order"]),
                 get_default_nested_dict(contig_read_level_purity, [c, "amb",
                                                                     "family"]),
                 get_default_nested_dict(contig_read_level_purity, [c, "amb",
                                                                     "genus"]),
                 get_default_nested_dict(contig_read_level_purity, [c, "amb",
                                                                     "species"]),
                 get_default_nested_dict(contig_read_level_purity, [c, "amb",
                                                                     "topname"]),
                 get_default_nested_dict(contig_read_level_purity, [c, "amb",
                                                                    "lca"],
                                         map_values=dict(topname="strain",
                                                         no_rank="life"))
                 ]) + "\n")

    # Make plots genome contig coverage
    a = numpy.genfromtxt(
        outdir + '/genome-contig-coverage.tsv', names=True, dtype=None)
    make_dir(outdir + '/plots')
    plt.plot(a['GC_content'], a['read_cov_ratio'], '.')
    plt.xlabel('GC content of genome')
    plt.ylabel('Read coverage ratio of genome')
    plt.savefig(outdir + '/plots/GC_content_vs_read_cov_ratio.pdf')
    plt.savefig(outdir + '/plots/GC_content_vs_read_cov_ratio.svg')
    plt.savefig(outdir + '/plots/GC_content_vs_read_cov_ratio.png')
    plt.clf()
    plt.plot(a['GC_content'], a['genome_contig_cov_ratio'], '.')
    plt.xlabel('GC content of genome')
    plt.ylabel('Contig coverage ratio of genome')
    plt.savefig(outdir + '/plots/GC_content_vs_genome_contig_cov_ratio.pdf')
    plt.savefig(outdir + '/plots/GC_content_vs_genome_contig_cov_ratio.svg')
    plt.savefig(outdir + '/plots/GC_content_vs_genome_contig_cov_ratio.png')
    plt.clf()
    plt.plot(a['read_cov_ratio'], a['genome_contig_cov_ratio'], '.')
    plt.xlabel('Read coverage ratio of genome')
    plt.ylabel('Genome contig coverage ratio of genome')
    plt.savefig(
        outdir + '/plots/read_cov_ratio_vs_genome_contig_cov_ratio.pdf')
    plt.savefig(
        outdir + '/plots/read_cov_ratio_vs_genome_contig_cov_ratio.svg')
    plt.savefig(
        outdir + '/plots/read_cov_ratio_vs_genome_contig_cov_ratio.png')
    plt.clf()

    # Make plots contig purity
    a = numpy.genfromtxt(outdir + '/contig-purity.tsv', names=True,
                         dtype=None, missing_values="-")
    plt.plot(a['contig_length'], a['unamb_read_level_purity'], '.')
    plt.xlabel('Contig length')
    plt.ylabel('Unambiguous read level purity of contig')
    plt.savefig(outdir + '/plots/contig_length_vs_unamb_read_level_purity.pdf')
    plt.savefig(outdir + '/plots/contig_length_vs_unamb_read_level_purity.svg')
    plt.savefig(outdir + '/plots/contig_length_vs_unamb_read_level_purity.png')
    plt.clf()
    plt.plot(a['contig_length'], a['amb_read_level_purity'], '.')
    plt.xlabel('Contig length')
    plt.ylabel('Ambiguous read level purity of contig')
    plt.savefig(outdir + '/plots/contig_length_vs_amb_read_level_purity.pdf')
    plt.savefig(outdir + '/plots/contig_length_vs_amb_read_level_purity.svg')
    plt.savefig(outdir + '/plots/contig_length_vs_amb_read_level_purity.png')
    plt.clf()
    plt.plot(a['contig_length'], a['max_aln_purity'], '.')
    plt.xlabel('Contig length')
    plt.ylabel('Alignment purity of contig')
    plt.savefig(outdir + '/plots/contig_length_vs_max_aln_purity.pdf')
    plt.savefig(outdir + '/plots/contig_length_vs_max_aln_purity.svg')
    plt.savefig(outdir + '/plots/contig_length_vs_max_aln_purity.png')
    plt.clf()
    plt.plot(a['unamb_read_level_purity'], a['max_aln_purity'], '.')
    plt.xlabel('Unambiguous read level purity of contig')
    plt.ylabel('Alignment purity of contig')
    plt.savefig(outdir + '/plots/unamb_read_level_purity_vs_max_aln_purity.pdf')
    plt.savefig(outdir + '/plots/unamb_read_level_purity_vs_max_aln_purity.svg')
    plt.savefig(outdir + '/plots/unamb_read_level_purity_vs_max_aln_purity.png')
    plt.clf()
    plt.plot(a['amb_read_level_purity'], a['max_aln_purity'], '.')
    plt.xlabel('Ambiguous read level purity of contig')
    plt.ylabel('Alignment purity of contig')
    plt.savefig(outdir + '/plots/amb_read_level_purity_vs_max_aln_purity.pdf')
    plt.savefig(outdir + '/plots/amb_read_level_purity_vs_max_aln_purity.svg')
    plt.savefig(outdir + '/plots/amb_read_level_purity_vs_max_aln_purity.png')
    plt.clf()

    chimer = a['max_aln_purity'] < 0.95
    d = defaultdict(int)
    for lca in a[chimer]['amb_lca']:
        if lca != "-":
            d[lca] += 1
    TAX_ORDER = ["strain", "species", "genus", "family", "order", "class",
                 "phylum", "superphylum", "superkingdom", "life"]
    bar_plot(d, TAX_ORDER, outdir=outdir)

    # Output plots in HTML
    sdir = os.path.dirname(os.path.realpath(__file__))
    template = open(sdir + '/validate-template.html').read()
    with open(outdir + '/index.html', 'w') as fh:
        fh.write(template.format(asmtsv='asm-stats.tsv',
                                 gcctsv='genome-contig-coverage.tsv',
                                 cptsv='contig-purity.tsv',
                                 plotgc1='plots/GC_content_vs_read_cov_ratio.png',
                                 plotgc2='plots/GC_content_vs_genome_contig_cov_ratio.png',
                                 plotgc3='plots/read_cov_ratio_vs_genome_contig_cov_ratio.png',
                                 plotcp1='plots/contig_length_vs_unamb_read_level_purity.png',
                                 plotcp2='plots/contig_length_vs_amb_read_level_purity.png',
                                 plotcp3='plots/contig_length_vs_max_aln_purity.png',
                                 plotcp4='plots/unamb_read_level_purity_vs_max_aln_purity.png',
                                 plotcp5='plots/amb_read_level_purity_vs_max_aln_purity.png',
                                 ))


def bar_plot(d, bars, outdir='.'):
    N = len(bars)
    ind = numpy.arange(N)
    width = 0.35
    plt.title("LCA of reads for contigs with alignment purity < 0.95", size='x-small')

    # add bars
    sum_values = sum(d.values())
    plt.bar(ind, [d[k] * 100 / sum_values for k in bars], width, color='red')
    # Final bar is not really visible, increase xlim
    xmin, xmax = plt.xlim()
    plt.xlim(xmin=xmin-1)
    plt.xlim(xmax=xmax+1)

    # axis setup
    labels = []
    for k in bars:
        labels.append('%s\n%d' % (k, d[k]))
    plt.xticks(ind + width / 2., bars, size='xx-small', rotation=17)
    tick_range = numpy.arange(0, 110, 10)
    plt.yticks(tick_range, size='xx-small')
    formatter = plt.FixedFormatter([str(x) for x in tick_range])
    plt.gca().yaxis.set_major_formatter(formatter)
    plt.gca().yaxis.grid(which='major')

    plt.savefig(outdir + "/plots/barplot.png")
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("coords", help="Output from nucmer show-coords")
    parser.add_argument("refstats", help="Reference stats tsv file\n")
    parser.add_argument("refphyl", help="Reference Phylogeny tsv file\n")
    parser.add_argument(
        "contigs", help="Contigs of the assembly in fasta format\n")
    parser.add_argument(
        "bamref", help="BAM file of the reads mapped against the reference\n")
    parser.add_argument(
        "bamasm", help="BAM file of the reads mapped against the contigs\n")
    parser.add_argument("outdir", help="Output directory\n")
    args = parser.parse_args()
    main(args.coords,
         args.refstats,
         args.refphyl,
         args.contigs,
         args.bamref,
         args.bamasm,
         args.outdir)
