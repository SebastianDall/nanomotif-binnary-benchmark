import os
import glob
import datetime
import polars as pl

SNAKEDIR = os.path.dirname(workflow.snakefile)
configfile: os.path.join(SNAKEDIR, "..", "config", "config-include.yaml")


# Extract sample names from the YAML configuration
# BENCHMARKS = [benchmark for benchmark in config['benchmarks'].keys()]
DATETIME = datetime.datetime.now().strftime('%Y%m%d_%H%M')
BASELINE = config["baseline_dir"]
OUTDIR = config["outdir"]

include: "rules/motif_discovery.smk"


rule all:
    input:
        expand(
            os.path.join(OUTDIR, "{sample}", "include.done"),
            sample=config["samples"].keys()
        ),
        expand(
            os.path.join(BASELINE, "include_files", "{sample}", f"motifs-scored-read-methylation-{config['min_valid_read_cov']}.tsv"),
            sample = config["samples"].keys()
        )



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

rule create_motifs_scored_file:
    input:
        a = lambda wildcards: config["samples"][wildcards.sample]["assembly"],
        p = lambda wildcards: config["samples"][wildcards.sample]["pileup"],
        b = lambda wildcards: config["samples"][wildcards.sample]["bin_motifs"],
    output:
        m = os.path.join(BASELINE, "include_files", "{sample}", f"motifs-scored-read-methylation-{config['min_valid_read_cov']}.tsv")
    threads:
        50
    resources:
        mem = "200G",
        walltime = "08:00:00",
        nodetype = "high-mem",
    params:
        MU = os.path.join(SNAKEDIR, "src", "methylation_utils")
    run:
        from pymethylation_utils.utils import run_methylation_utils
        import sys

        motifs = pl.read_csv(input[2], separator = "\t")
        motifs = motifs\
            .with_columns((pl.col("motif") + "_" + pl.col("mod_type") + "_" + pl.col("mod_position").cast(pl.Utf8)).alias("motif_mod"))\
            .get_column("motif_mod").unique()

        print("Running methylation_utils")
        return_code = run_methylation_utils(
            pileup = input[1],
            assembly = input[0],
            motifs = motifs,
            threads = threads,
            min_valid_read_coverage = config["min_valid_read_cov"],
            output = output[0]
        )

        if return_code != 0:
            sys.exit(1)

rule symlink_motifs_scored_file:
    input:
        m = os.path.join(BASELINE, "include_files", "{sample}", f"motifs-scored-read-methylation-{config['min_valid_read_cov']}.tsv")
    output:
        m = touch(os.path.join(OUTDIR,"{sample}","{benchmark}", "motifs-read-methylation-copied.done"))
    threads:
        1
    resources:
        mem = "4G",
        walltime = "00:02:00",
        nodetype = config["partition"],
    shell:
        """
            ln -s "$(realpath {input.m})" {OUTDIR}/{wildcards.sample}/{wildcards.benchmark}/motifs-scored-read-methylation.tsv
        """

checkpoint create_benchmark_contig_bin_files:
    conda: "nanomotif_dev"
    input:
        m = os.path.join(BASELINE, "include_files", "{sample}", f"motifs-scored-read-methylation-{config['min_valid_read_cov']}.tsv"),
        c = lambda wildcards: config["samples"][wildcards.sample]["contig_bin"]
    output:
        directory(os.path.join(BASELINE, "benchmark_files", "{sample}", "include"))
    threads:
        1
    resources:
        mem = "10G",
        walltime = "00:15:00",
        nodetype = config["partition"],
    params:
        PY = os.path.join(SNAKEDIR, "src", "shuffle_contigs.py"),
        n_benchmarks = config["n_benchmarks"]
    shell:
        """
            python {params.PY} {input.c} {input.m} {params.n_benchmarks} include {BASELINE}/benchmark_files/{wildcards.sample}/
        """

        

    
        
rule symlink_contig_bin_dev:
    input:
        os.path.join(BASELINE, "benchmark_files", "{sample}", "include", "{benchmark}.tsv")
    output:
        os.path.join(OUTDIR,"{sample}","{benchmark}",  "contig_bin.tsv"),
    threads:
        1
    resources:
        mem = "1G",
        walltime = "2:00",
        nodetype = config["partition"],
    shell:
        """
        ln -s "$(realpath {input})" {output}
        """
        

rule include_dev:
    conda:
        "nanomotif_dev"
    input:
        a = lambda wildcards: config["samples"][wildcards.sample]["assembly"],
        p = lambda wildcards: config["samples"][wildcards.sample]["pileup"],
        b = lambda wildcards: config["samples"][wildcards.sample]["bin_motifs"],
        bin_con = os.path.join(BASELINE, "include_files", "empty_bin_contamination.tsv"),
        c = os.path.join(OUTDIR,"{sample}","{benchmark}",  "contig_bin.tsv"),
        m = os.path.join(OUTDIR,"{sample}","{benchmark}", "motifs-read-methylation-copied.done")
    output:
        os.path.join(OUTDIR,"{sample}","{benchmark}",  "include_contigs.tsv"),
    threads:
        5
    resources:
        mem = "10G",
        walltime = "08:00:00",
        nodetype = config["partition"],
    shell:
        """
        nanomotif include_contigs \
            --contamination_file {input.bin_con} \
            --threads {threads} \
            --pileup {input.p} \
            --assembly {input.a} \
            --bin_motifs {input.b} \
            --contig_bins {input.c} \
            --out {OUTDIR}/{wildcards.sample}/{wildcards.benchmark}/
        """


def find_benchmarks(wildcards):
    cp = checkpoints.create_benchmark_contig_bin_files.get(**wildcards).output[0]

    return expand(
        os.path.join(OUTDIR, "{sample}", "{benchmark}", "include_contigs.tsv"),
        sample = wildcards.sample, benchmark = glob_wildcards(os.path.join(cp, "{benchmark}.tsv")).benchmark
    )

rule aggregate_benchmarks:
    input:
        find_benchmarks
    output:
        touch(os.path.join(OUTDIR, "{sample}", "include.done"))
    threads:
        1
    resources:
        mem = "2G",
        walltime = "00:02:00",
        nodetype = "default-op",
