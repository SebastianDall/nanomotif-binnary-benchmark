import os
import polars as pl
import pandas as pd
import random
import copy

# PERCENT_OF_BINS_WITH_CONTAMINATION = [0.1, 0.3, 0.5, 0.7, 0.9]
# NUMBER_OF_CONTAMINANTS = [1, 2, 3]

def associate_bin_with_contigs(contig_bin):
    bin_contigs = {}
    for bin in contig_bin["bin"].unique():
        bin_contigs[bin] = contig_bin[contig_bin["bin"] == bin]["contig"].to_list()

    return bin_contigs
    


def shuffle_contigs_between_bins(bin_contigs, num_contigs_to_shuffle):
    contigs_to_move = {}

    for bin in bin_contigs.keys():
        contigs = bin_contigs[bin]
        # INFO: Since we are no longer removing contigs. Bigger contigs can also be copied to other bins.
        # if len(contigs) <= num_contigs_to_shuffle:
        #     continue

        print("Processing bin:", bin)
        print("Processing contigs:", contigs)
        contigs_for_shuffling = random.sample(contigs, num_contigs_to_shuffle)

        contigs_to_move[bin] = contigs_for_shuffling

        # INFO: Removing contigs at random can create some weird cases where the biggest contig in a bin is moved. This will mean a bin could be a very low quality bin that we are trying to assess.
        # # remove contigs to be shuffled.
        # bin_contigs[bin] = [c for c in contigs if c not in contigs_for_shuffling]

    shuffled_items = list(contigs_to_move.items())
    random.shuffle(shuffled_items)
    
    bin_wo_contaminants = list(bin_contigs.keys())
    for source_bin, contigs in shuffled_items:
        for contig in contigs:
            possible_bins_for_assignment = [b for b in contigs_to_move.keys() if b != source_bin and "coli" not in b and b in bin_wo_contaminants]
            if len(possible_bins_for_assignment) == 0:
                print(f"No possible bin for assigning {contig} from {source_bin}")
                break
            
            new_bin = random.choice(possible_bins_for_assignment)
            bin_contigs[new_bin].append("contaminant_" + contig)
            bin_wo_contaminants.remove(new_bin)

    return bin_contigs

# def create_pileup(pileup, contaminated_bin_contigs):
    

def remove_contigs_from_bin(bin_contigs, num_contigs_to_remove):
    contigs_to_remove = {}

    for bin in bin_contigs.keys():
        contigs = bin_contigs[bin]
        if len(contigs) <= num_contigs_to_remove:
            continue

        contigs_for_removal = random.sample(contigs, num_contigs_to_remove)

        contigs_to_remove[bin] = contigs_for_removal

        # remove contigs to be shuffled.
        bin_contigs[bin] = [c for c in contigs if c not in contigs_for_removal]

    return bin_contigs

def write_contig_bin_file(bin_contigs, filename):
    with open(filename, "w") as f:
        for bin, contigs in bin_contigs.items():
            for contig in contigs:
                f.write(f"{contig}\t{bin}\n")    
    print(f"Created file {filename}")


if __name__ == '__main__':
    random.seed(42)
    contig_bin = pd.read_csv("data/datasets/simulated_4_lognormal_20-100x/contig_mapping/mapped_contig_bin.tsv", sep = "\t", names=["contig", "bin"])

    contig_methylation = pl.read_csv("data/datasets/simulated_4_lognormal_20-100x/motifs-scored-read-methylation.tsv", separator = "\t")

    contigs_w_methylation = contig_methylation\
        .filter(pl.col("N_motif_obs") > 8)\
        .get_column("contig").unique()

    contig_bin = contig_bin[contig_bin["contig"].isin(contigs_w_methylation)]
    contigs_in_bin = associate_bin_with_contigs(contig_bin)


    for i in range(0,10):
        print(f"Shuffling benchmark {i}")
        shuffled_contigs = shuffle_contigs_between_bins(copy.deepcopy(contigs_in_bin), 1)
        write_contig_bin_file(shuffled_contigs, f"files/benchmarks/simulated_4_lognormal_20-100x/benchmark_{i}_shuffle_1_contig_w_signal.tsv")
        
        shuffled_contigs = remove_contigs_from_bin(copy.deepcopy(contigs_in_bin), 1)
        write_contig_bin_file(shuffled_contigs, f"files/benchmarks/simulated_4_lognormal_20-100x/benchmark_{i}_remove_1_contig_w_signal.tsv")
        
