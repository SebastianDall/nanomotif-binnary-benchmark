import os
import glob
import datetime
import polars as pl

SNAKEDIR = os.path.dirname(workflow.snakefile)
configfile: os.path.join(SNAKEDIR, "..", "config", "config-motif_discovery.yaml")


# Extract sample names from the YAML configuration
# BENCHMARKS = [benchmark for benchmark in config['benchmarks'].keys()]
DATETIME = datetime.datetime.now().strftime('%Y%m%d_%H%M')
OUTDIR = config["outdir"]

rule all:
    input:
        expand(
        os.path.join(OUTDIR, "{sample}", "nanomotif_0.4.16", "bin-motifs.tsv"),
            sample=config["samples"].keys()
        ),
        

rule nanomotif_motif_discovery:
    conda:
        "envs/nanomotif.yaml"
    input:
        a = lambda wildcards: config["samples"][wildcards.sample]["assembly"],
        p = lambda wildcards: config["samples"][wildcards.sample]["pileup"],
        bins = lambda wildcards: config["samples"][wildcards.sample]["contig_bin"],
    output:
        o = os.path.join(OUTDIR, "{sample}", "nanomotif_0.4.16", "bin-motifs.tsv"),
        motifs = os.path.join(OUTDIR, "{sample}", "nanomotif_0.4.16","motifs.tsv"),
        m = os.path.join(OUTDIR, "{sample}", "nanomotif_0.4.16","motifs-scored.tsv"),
    threads:
        40
    resources:
        mem = "80G",
        walltime = "1-00:00:00",
        nodetype = "high-mem",
    shell:
        """
        nanomotif motif_discovery \
            -t {threads} \
            --out {OUTDIR}/{wildcards.sample}/nanomotif_0.4.16/ \
            {input.a} {input.p} {input.bins}
        """

