import os
import polars as pl
import pandas as pd
import random

PERCENT_OF_BINS_WITH_CONTAMINATION = [0.1, 0.3, 0.5, 0.7, 0.9]
NUMBER_OF_CONTAMINANTS = [1, 2, 3]

contig_bin = pd.read_csv("data/datasets/simulated_3_lognormal/contig_mapping/mapped_contig_bin.tsv", sep = "\t")
bin_contigs = {}
for bin in contig_bin["bin"].unique():
    bin_contigs[bin] = contig_bin[contig_bin["bin"] == bin]["contig"].to_list()


def shuffle_contigs_between_bins(bin_contigs, num_contigs_to_shuffle):
    contigs_to_move = {}

    for bin in bin_contigs.keys():
        contigs = bin_contigs[bin]
        if len(contigs) <= num_contigs_to_shuffle:
            continue

        contigs_for_shuffling = random.sample(contigs, num_contigs_to_shuffle)

        contigs_to_move[bin] = contigs_for_shuffling

        # remove contigs to be shuffled.
        bin_contigs[bin] = [c for c in contigs if c not in contigs_for_shuffling]

    for source_bin, contigs in contigs_to_move.items():
        for contig in contigs:
            possible_bins_for_assignment = [b for b in contigs_to_move.keys() if b != source_bin]
            new_bin = random.choice(possible_bins_for_assignment)
            bin_contigs[new_bin].append(contig)

    return bin_contigs



def write_contig_bin_file(bin_contigs, filename):
    with open(filename, "w") as f:
        for bin, contigs in bin_contigs.items():
            for contig in contigs:
                f.write(f"{contig}\t{bin}\n")    


new_contig_bin = shuffle_contigs_between_bins(bin_contigs, 1)
write_contig_bin_file(new_contig_bin, "test.tsv")
