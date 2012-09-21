#!/usr/bin/awk -f
# Parses a coords file from nucmer and outputs purity information of the entire
# assembly.
BEGIN {
    if (!name) name="-"
    if (!kmer_type) kmer_type="-"
    if (!asm_type) asm_type="-"
    if (!kmer_size) kmer_size="-"
    if (!kmin) kmin="-"
    if (!kmax) kmax="-"

    cut_off = 100

    # L50 variables
    sum_bases = 0
    nr_contigs = 0
}
{
    # First determine fasta file lengths
    if (NR == FNR) { 
        if ($0 !~ /^>/ && length($0) >= cut_off) {
            # Save the contig lengths in a way that allows awk to sort them
            all_contig_lengths[nr_contigs++] = sprintf("%12s", length($0))
            sum_bases += length($0)
        }
        next
    }
    # Store purity and length of each contig
    if (strtonum($9) >= cut_off && strtonum($11) > max_qcov[$13]) {
        max_qcov[$13] = strtonum($11)
        contig_length[$13] = strtonum($9)
    }
    #
    if (!($12 in ref_lengths)) {
        ref_lengths[$12] = strtonum($8)
    }
}
END {
    # Calculate L50 
    asort(all_contig_lengths)
    cumsum = 0
    for (i = nr_contigs; i > 0; i--) {
        cumsum += strtonum(all_contig_lengths[i])
        if (cumsum >= sum_bases / 2) {
            n50 = nr_contigs + 1 - i
            l50 = strtonum(all_contig_lengths[i])
            break
        }
    }
    # Calculate purity
    sum_max_qcov_bases = 0
    nr_contigs_mapping = 0
    for (contig in max_qcov) {
        sum_max_qcov_bases += contig_length[contig] * max_qcov[contig] / 100
        nr_contigs_mapping++
    }
    # Sum reference lengths
    sum_ref_lengths = 0
    for (ref in ref_lengths) {
        sum_ref_lengths += ref_lengths[ref]
    }
    
    print "name","purity","l50","n50","trim_n","max_contig_length","cut_off","trim_n_mapping","sum_ref_lengths","sum_max_qcov_bases","ratio","sum_bases","asm_type","kmer_type","kmer_size","kmin","kmax"
    print name,sum_max_qcov_bases / sum_bases,l50,n50,nr_contigs,all_contig_lengths[nr_contigs],cut_off,nr_contigs_mapping,sum_ref_lengths,sum_max_qcov_bases,sum_max_qcov_bases/ sum_ref_lengths,sum_bases,asm_type,kmer_type,kmer_size,kmin,kmax

    # Print contig_length versus purity
    if (contig_out) {
        print "contig","contig_length","max_qcov" > contig_out
        for (contig in max_qcov) {
            print contig,contig_length[contig],max_qcov[contig] > contig_out
        }
    }
}
