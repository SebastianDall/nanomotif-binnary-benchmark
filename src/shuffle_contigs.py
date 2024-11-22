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
    


def shuffle_contigs_between_bins(contig_bin, contigs_to_shuffle, num_contigs_to_shuffle):
    new_contig_bin = copy.deepcopy(contig_bin)
    
    bins_wo_contaminants = list(contig_bin.keys())
    bins_wo_contaminants = [b for b in bins_wo_contaminants if "coli" not in b]

    bins_where_contaminants_have_been_taken = set()
    contigs_used_as_contaminants = set()
    
    while bins_wo_contaminants:
        print("Remaining bins without contaminants: ", len(bins_wo_contaminants))
        print(bins_wo_contaminants)

        bins_availble_for_taking_contamination = [b for b in contigs_to_shuffle.keys() if b not in bins_where_contaminants_have_been_taken]
        if not bins_availble_for_taking_contamination:
            bins_availble_for_taking_contamination = list(bins_where_contaminants_have_been_taken)

        bin_to_take_from = random.choice(bins_availble_for_taking_contamination)
        contigs_available = contigs_to_shuffle[bin_to_take_from]
        contigs_available = [c for c in contigs_available if c not in contigs_used_as_contaminants]
        if not contigs_available:
            bins_where_contaminants_have_been_taken.add(bin_to_take_from)
            continue

        contig_selected = random.choice(contigs_available)
        bin_to_assign_to = random.choice(bins_wo_contaminants)

        if not contig_selected in contig_bin[bin_to_assign_to]:
            new_contig_bin[bin_to_assign_to].append("contamination_" + contig_selected)
            bins_wo_contaminants.remove(bin_to_assign_to)
            contigs_used_as_contaminants.add(contig_selected)
            bins_where_contaminants_have_been_taken.add(bin_to_take_from)

    return new_contig_bin

def remove_contigs_from_bin(bin_contigs, contigs_available_for_removal, num_contigs_to_remove):
    contigs_to_remove = {}

    available_bins =contigs_available_for_removal.keys()
    for bin in available_bins:
        contigs = contigs_available_for_removal[bin]
        if not contigs:
            continue

        contigs_for_removal = random.choice(contigs)

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
    contig_bin = pl.read_csv("data/datasets/simulated_4_lognormal_20-100x/contig_mapping/mapped_contig_bin.tsv", 
                             separator = "\t", has_header=False, new_columns=["contig", "bin"])
    contig_methylation = pl.read_csv("data/datasets/simulated_4_lognormal_20-100x/motifs-scored-read-methylation.tsv", separator = "\t")

    contigs_w_methylation = contig_methylation\
        .filter(pl.col("N_motif_obs") > 8)

    bin_observations = contigs_w_methylation.join(contig_bin, how = "left", on = ["contig"])\
        .filter((pl.col("motif") == "GATC") & (pl.col("mod_type") == "a"))\
        .group_by(["bin", "motif"])\
        .agg(
            pl.col("N_motif_obs").sum().alias("N_bin_motif_obs")
        )
        
    smaller_contigs_w_methylation = contigs_w_methylation\
        .filter((pl.col("motif") == "GATC") & (pl.col("mod_type") == "a"))\
        .join(contig_bin, how = "left", on = ["contig"])\
        .join(bin_observations, how="left", on = ["bin"])\
        .with_columns(
            obs_percentage = pl.col("N_motif_obs") / pl.col("N_bin_motif_obs")
        )

    print(smaller_contigs_w_methylation.filter(pl.col("contig").is_in(["contig_329", "contig_63"])))
        
    smaller_contigs_w_methylation = smaller_contigs_w_methylation.filter(pl.col("obs_percentage") < 0.5).get_column("contig")
        

    contig_bin = contig_bin.to_pandas()
    contigs_to_shuffle = contig_bin[contig_bin["contig"].isin(smaller_contigs_w_methylation)]
    contigs_to_shuffle = associate_bin_with_contigs(contigs_to_shuffle)

    print(contigs_to_shuffle)
    contig_bin_dict = associate_bin_with_contigs(contig_bin)

    for i in range(0,10):
        print(f"Shuffling benchmark {i}")
        shuffled_contigs = shuffle_contigs_between_bins(copy.deepcopy(contig_bin_dict), copy.deepcopy(contigs_to_shuffle), 1)
        write_contig_bin_file(shuffled_contigs, f"files/benchmarks/simulated_4_lognormal_20-100x/benchmark_{i}_shuffle_1_contig_w_signal.tsv")
        
        shuffled_contigs = remove_contigs_from_bin(copy.deepcopy(contig_bin_dict), copy.deepcopy(contigs_to_shuffle), 1)
        write_contig_bin_file(shuffled_contigs, f"files/benchmarks/simulated_4_lognormal_20-100x/benchmark_{i}_remove_1_contig_w_signal.tsv")
        
