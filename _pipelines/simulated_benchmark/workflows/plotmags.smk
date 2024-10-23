import os
import glob
import datetime

SNAKEDIR = os.path.dirname(workflow.snakefile)
configfile: os.path.join(SNAKEDIR, "..", "config", "config-contamination.yaml")

BENCHMARK_DIR = config["outdir"]
BASELINE_DIR = config["baseline_dir"]
MODE = config["mode"]
SAMPLES = config["samples"].keys()

BASELINE_DIR_MODE = os.path.join(BASELINE_DIR, MODE)

OUTDIRS = [BENCHMARK_DIR, BASELINE_DIR_MODE]

def benchmark_inputs():
    inputs = []
    for dir in OUTDIRS:
        for sample in SAMPLES:
            for benchmark in config["samples"][sample]["benchmarks"].keys():
                inputs.append(os.path.join(dir, sample, benchmark, "MAGs", "done.txt"))
    return inputs

rule all:
    input:
        benchmark_inputs()


nanomotif_config = {}
nanomotif_config[BENCHMARK_DIR] = {
    "mean_methylation_cutoff": f'--mean_methylation_cutoff {config["mean_methylation_cutoff"]}',
    "n_motif_contig_cutoff": f'--n_motif_contig_cutoff {config["n_motif_contig_cutoff"]}',
    "n_motif_bin_cutoff": f'--n_motif_bin_cutoff {config["n_motif_bin_cutoff"]}'    
}
nanomotif_config[BASELINE_DIR_MODE] = {
    "mean_methylation_cutoff":  f'--mean_methylation_cutoff {0.25}',
    "n_motif_contig_cutoff":  f'--n_motif_contig_cutoff {10}',
    "n_motif_bin_cutoff": f'--n_motif_bin_cutoff {500}'
}

rule plot_MAGs:
    input:
        motifs_scored = os.path.join(config["baseline_dir"], "motif_discovery", "{sample}","original_contig_bin", "motifs-scored.tsv"),
        bin_motifs = os.path.join(config["baseline_dir"], "motif_discovery", "{sample}","original_contig_bin", "bin-motifs.tsv"),
        contig_bins = lambda wildcards: config["samples"][wildcards.sample]["benchmarks"][wildcards.benchmark],
        contig_bins_truth = lambda wildcards: config["samples"][wildcards.sample]["benchmarks"]["original_contig_bin"],
        bin_contamination = os.path.join("{dir}", "{sample}", "{benchmark}","bin_contamination.tsv"),
    output:
        touch(os.path.join("{dir}",  "{sample}", "{benchmark}", "MAGs", "done.txt"))
    threads:
        1
    resources:
        mem = "2G",
        walltime = "00:25:00",
        nodetype = config["partition"],
    params:
        script = os.path.join(SNAKEDIR, "src","plot_MAGs.R"),
        mean_methylation_cutoff = lambda wildcards: nanomotif_config[wildcards.dir]["mean_methylation_cutoff"], 
        n_motif_contig_cutoff = lambda wildcards: nanomotif_config[wildcards.dir]["n_motif_contig_cutoff"], 
        n_motif_bin_cutoff = lambda wildcards: nanomotif_config[wildcards.dir]["n_motif_bin_cutoff"],
    shell:
        """
        Rscript {params.script} --motifs_scored {input.motifs_scored} \
            --bin_motifs {input.bin_motifs} \
            --contig_bins {input.contig_bins} \
            --contig_bins_truth {input.contig_bins_truth} \
            --bin_contamination {input.bin_contamination} \
            {params.mean_methylation_cutoff} \
            {params.n_motif_contig_cutoff} \
            {params.n_motif_bin_cutoff} \
            --output {wildcards.dir}/{wildcards.sample}/{wildcards.benchmark}/MAGs  
        """
