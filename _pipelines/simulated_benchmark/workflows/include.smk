import os
import glob
import datetime

SNAKEDIR = os.path.dirname(workflow.snakefile)
configfile: os.path.join(SNAKEDIR, "..", "config", "config-include.yaml")


# Extract sample names from the YAML configuration
# BENCHMARKS = [benchmark for benchmark in config['benchmarks'].keys()]
DATETIME = datetime.datetime.now().strftime('%Y%m%d_%H%M')
BASELINE = config["baseline_dir"]
OUTDIR = config["outdir"]

include: "rules/motif_discovery.smk"

def benchmark_inputs():
    inputs = []

    for sample in config["samples"].keys():
        for benchmark in config["samples"][sample]["benchmarks"].keys():
            inputs.append(os.path.join(BASELINE, "include", sample, benchmark, "include_contigs.tsv"))
            inputs.append(os.path.join(OUTDIR, sample, benchmark, "include_contigs.tsv"))
    return inputs


rule all:
    input:
        benchmark_inputs(),
        expand(os.path.join(BASELINE, "motif_discovery", "{sample}",  "original_contig_bin","motifs-scored.tsv"),
               sample = config["samples"].keys()),
       


rule copy_contig_bin_baseline:
    input:
        lambda wildcards: config["samples"][wildcards.sample]["benchmarks"][wildcards.benchmark],
    output:
        os.path.join(BASELINE, "include", "{sample}", "{benchmark}", "contig_bin.tsv"),
    threads:
        1
    resources:
        mem = "1G",
        walltime = "2:00",
        nodetype = config["partition"],
    shell:
        "cp {input} {output}"

rule copy_contig_bin_dev:
    input:
        lambda wildcards: config["samples"][wildcards.sample]["benchmarks"][wildcards.benchmark],
    output:
        os.path.join(OUTDIR,"{sample}","{benchmark}",  "contig_bin.tsv"),
    threads:
        1
    resources:
        mem = "1G",
        walltime = "2:00",
        nodetype = config["partition"],
    shell:
        "cp {input} {output}"
        

rule create_empty_contamination_file:
    output:
        os.path.join(BASELINE, "include_files", "empty_bin_contamination.tsv")
    threads: 1
    resources:
        mem = "2G",
        walltime = "2:00",
        nodetype = config["partition"]
    shell:
        """
        echo -e "bin\tbin_contig_compare\tbinary_methylation_missmatch_score\tnon_na_comparisons\tcontig" > {output}
        """

rule nanomotif_include:
    conda:
        "envs/nanomotif.yaml"
    input:
        m = os.path.join(BASELINE, "motif_discovery", "{sample}", "original_contig_bin", "motifs-scored.tsv"),
        b = os.path.join(BASELINE, "motif_discovery", "{sample}", "original_contig_bin", "bin-motifs.tsv"),
        c = os.path.join(BASELINE, "include", "{sample}", "{benchmark}", "contig_bin.tsv"),
        bin_con = os.path.join(BASELINE, "include_files", "empty_bin_contamination.tsv")
    output:
        os.path.join(BASELINE, "include", "{sample}", "{benchmark}", "include_contigs.tsv"),
    threads:
        20
    resources:
        mem = "40G",
        walltime = "08:00:00",
        nodetype = config["partition"],
    params:
    shell:
        """
        nanomotif include_contigs \
            --contamination_file {input.bin_con} \
            --motifs_scored {input.m} \
            --bin_motifs {input.b} \
            --contig_bins {input.c} \
            --out {BASELINE}/detect_contamination/{wildcards.sample}/{wildcards.benchmark}         
        """

rule nanomotif_include_dev:
    conda:
        "nanomotif_dev"
    input:
        a = lambda wildcards: config["samples"][wildcards.sample]["assembly"],
        p = lambda wildcards: config["samples"][wildcards.sample]["pileup"],
        b = os.path.join(BASELINE, "motif_discovery", "{sample}", "original_contig_bin", "bin-motifs.tsv"),
        bin_con = os.path.join(BASELINE, "include_files", "empty_bin_contamination.tsv"),
        c = os.path.join(OUTDIR,"{sample}","{benchmark}",  "contig_bin.tsv"),
    output:
        os.path.join(OUTDIR,"{sample}","{benchmark}",  "include_contigs.tsv"),
    threads:
        40
    resources:
        mem = "80G",
        walltime = "08:00:00",
        nodetype = config["partition"],
    params:
        mean_methylation_cutoff = config["mean_methylation_cutoff"],
        n_motif_contig_cutoff = config["n_motif_contig_cutoff"],
        n_motif_bin_cutoff = config["n_motif_bin_cutoff"],
    shell:
        """
        nanomotif include_contigs \
            --contamination_file {input.bin_con} \
            --threads {threads} \
            --pileup {input.p} \
            --assembly {input.a} \
            --bin_motifs {input.b} \
            --contig_bins {input.c} \
            --mean_methylation_cutoff {params.mean_methylation_cutoff} \
            --n_motif_contig_cutoff {params.n_motif_contig_cutoff} \
            --n_motif_bin_cutoff {params.n_motif_bin_cutoff} \
            --out {OUTDIR}/{wildcards.sample}/{wildcards.benchmark}/
        """
