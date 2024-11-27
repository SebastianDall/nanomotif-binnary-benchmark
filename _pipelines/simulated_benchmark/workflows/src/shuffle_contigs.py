import os
import polars as pl
import pandas as pd
import random
import copy
import sys


def associate_bin_with_contigs(contig_bin):
    bin_contigs = {}
    for bin in contig_bin["bin"].unique():
        bin_contigs[bin] = contig_bin[contig_bin["bin"] == bin]["contig"].to_list()

    return bin_contigs
    


def shuffle_contigs_between_bins(contig_bin, contigs_to_shuffle):
    new_contig_bin = copy.deepcopy(contig_bin)
    
    print(contig_bin.keys())
    bins_wo_contaminants = list(contig_bin.keys())

    bins_where_contaminants_have_been_taken = set()
    contigs_used_as_contaminants = set()
    
    i=0
    n=1
    while bins_wo_contaminants:
        if i == 50:
            sys.exit(1)
        print("Remaining bins without contaminants: ", len(bins_wo_contaminants))
        if len(bins_wo_contaminants) < 5:
            print(bins_wo_contaminants)
        print("bin where contamination has been taken")
        print(bins_where_contaminants_have_been_taken)

        bins_availble_for_taking_contamination = [b for b in contigs_to_shuffle.keys() if b not in bins_where_contaminants_have_been_taken]
        print("bins == bins available: ",bins_wo_contaminants == bins_availble_for_taking_contamination)
        print("n = ", n)
        if not bins_availble_for_taking_contamination or (bins_wo_contaminants == bins_availble_for_taking_contamination and n !=1 ):
            print("all bins have had a contig taken as contaminant")
            bins_availble_for_taking_contamination = list(bins_where_contaminants_have_been_taken)

        if len(bins_wo_contaminants) < 5:
            print("Bins availble for taking contamination")
            print(bins_availble_for_taking_contamination)
        

        bin_to_take_from = random.choice(bins_availble_for_taking_contamination)
        contigs_available = contigs_to_shuffle[bin_to_take_from]
        contigs_available = [c for c in contigs_available if c not in contigs_used_as_contaminants]
        if not contigs_available:
            bins_where_contaminants_have_been_taken.add(bin_to_take_from)
            continue

        
        print("bins without contaminants: ", bins_wo_contaminants)
        contig_selected = random.choice(contigs_available)
        bin_to_assign_to = random.choice(bins_wo_contaminants)

        assert contig_selected is not None, "contig selected is none"
        assert bin_to_assign_to is not None, "Bin is None"
        print("contig selected: ",contig_selected)
        print("Bin to assign: ", bin_to_assign_to)
        if not contig_selected in contig_bin[bin_to_assign_to]:
            new_contig_bin[bin_to_assign_to].append("contamination_" + contig_selected)
            new_contig_bin[bin_to_take_from].remove(contig_selected)
            bins_wo_contaminants.remove(bin_to_assign_to)
            contigs_used_as_contaminants.add(contig_selected)
            bins_where_contaminants_have_been_taken.add(bin_to_take_from)
            i = 0
            n = 2
        i+=1

    return new_contig_bin

def remove_contigs_from_bin(bin_contigs, contigs_available_for_removal):
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
    contig_bin_file = sys.argv[1]
    contig_bin = pl.read_csv(contig_bin_file, 
                             separator = "\t", has_header=False, new_columns=["contig", "bin"])

    motifs_scored_file = sys.argv[2]
    contig_methylation = pl.read_csv(motifs_scored_file, separator = "\t")

    contigs = contig_bin.get_column("contig")
    contigs_w_methylation = contig_methylation\
        .filter(pl.col("N_motif_obs") > 8)\
        .filter(pl.col("contig").is_in(contigs))

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

       
    smaller_contigs_w_methylation = smaller_contigs_w_methylation.filter(pl.col("obs_percentage") < 0.5).get_column("contig")
        

    contig_bin = contig_bin.to_pandas()
    contigs_to_shuffle = contig_bin[contig_bin["contig"].isin(smaller_contigs_w_methylation)]
    contigs_to_shuffle = associate_bin_with_contigs(contigs_to_shuffle)

    contig_bin_dict = associate_bin_with_contigs(contig_bin)

    n_benchmarks = int(sys.argv[3])
    command = sys.argv[4]
    outdir = sys.argv[5]
    outdir = os.path.join(outdir, command)
    os.makedirs(outdir, exist_ok=True)

    for i in range(0,n_benchmarks):
        if command == "contamination":
            print(f"Shuffling benchmark {i}")
            shuffled_contigs = shuffle_contigs_between_bins(copy.deepcopy(contig_bin_dict), copy.deepcopy(contigs_to_shuffle))
            write_contig_bin_file(shuffled_contigs, os.path.join(outdir, f"benchmark_{i}_shuffle_1_contig_w_signal.tsv"))
        elif command == "include":
            shuffled_contigs = remove_contigs_from_bin(copy.deepcopy(contig_bin_dict), copy.deepcopy(contigs_to_shuffle))
            write_contig_bin_file(shuffled_contigs, os.path.join(outdir, f"benchmark_{i}_remove_1_contig_w_signal.tsv"))
        else:
            print("Wrong command")
            sys.exit(1)
        
