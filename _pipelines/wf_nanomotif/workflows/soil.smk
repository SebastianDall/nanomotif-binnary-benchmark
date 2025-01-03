import os
import glob
import datetime
import polars as pl

SNAKEDIR = os.path.dirname(workflow.snakefile)
configfile: os.path.join(SNAKEDIR, "..", "config", "config-soil.yaml")


# Extract sample names from the YAML configuration
# BENCHMARKS = [benchmark for benchmark in config['benchmarks'].keys()]
DATETIME = datetime.datetime.now().strftime('%Y%m%d_%H%M')
OUTDIR = config["outdir"]

rule all:
    input:
        expand(
            os.path.join(OUTDIR, "{sample}", "nanomotif_binnary", f"motifs-scored-read-methylation-{config['min_valid_read_cov']}.tsv"),
            sample = config["samples"].keys()
        ),
        expand(
            os.path.join(OUTDIR,"{sample}","nanomotif_binnary",  "bin_contamination.tsv"),
            sample = config["samples"].keys()
        ),
        expand(
            os.path.join(OUTDIR,"{sample}","nanomotif_binnary",  "include_contigs.tsv"),
            sample = config["samples"].keys()
        ),
        expand(
            os.path.join(OUTDIR, "{sample}", "nanomotif_binnary", "checkm2_decon", "quality_report.tsv"),
            sample = config["samples"].keys()
        )

def threads(num):
    threads_from_gb = int(num.replace("G","")) // 10
    return min(max(40, threads_from_gb), 90)        


rule create_motifs_scored_file:
    input:
        a = lambda wildcards: config["samples"][wildcards.sample]["assembly"],
        p = lambda wildcards: config["samples"][wildcards.sample]["pileup"],
        b = lambda wildcards: config["samples"][wildcards.sample]["bin_motifs"],
    output:
        m = os.path.join(OUTDIR, "{sample}", "nanomotif_binnary", f"motifs-scored-read-methylation-{config['min_valid_read_cov']}.tsv"),
    threads:
        lambda wildcards: threads(config["samples"][wildcards.sample]["mem"])
    resources:
        mem = lambda wildcards: config["samples"][wildcards.sample]["mem"],
        walltime = "1-10:00:00",
        nodetype = config["highmem"],
    params:
        batches= lambda wildcards: config["samples"][wildcards.sample]["batches"]
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
            batches = params[0],
            output = output[0]
        )

        if return_code != 0:
            sys.exit(1)

rule symlink_motifs_scored_file:
    input:
        m = os.path.join(OUTDIR, "{sample}", "nanomotif_binnary", f"motifs-scored-read-methylation-{config['min_valid_read_cov']}.tsv"),
    output:
        m = touch(os.path.join(OUTDIR,"{sample}","nanomotif_binnary", "motifs-read-methylation-copied.done"))
    threads:
        1
    resources:
        mem = "4G",
        walltime = "00:02:00",
        nodetype = config["partition"],
    shell:
        """
            ln -sf "$(realpath {input.m})" {OUTDIR}/{wildcards.sample}/nanomotif_binnary/motifs-scored-read-methylation.tsv
        """

rule nanomotif_contamination:
    conda:
        "nanomotif_dev"
    input:
        a = lambda wildcards: config["samples"][wildcards.sample]["assembly"],
        p = lambda wildcards: config["samples"][wildcards.sample]["pileup"],
        b = lambda wildcards: config["samples"][wildcards.sample]["bin_motifs"],
        c = lambda wildcards: config["samples"][wildcards.sample]["contig_bin"],
        m = os.path.join(OUTDIR,"{sample}","nanomotif_binnary", "motifs-read-methylation-copied.done")
    output:
        os.path.join(OUTDIR,"{sample}","nanomotif_binnary",  "bin_contamination.tsv"),
    threads:
        10
    resources:
        mem = "60G",
        walltime = "08:00:00",
        nodetype = config["partition"],
    shell:
        """
        nanomotif detect_contamination \
            --threads {threads} \
            --pileup {input.p} \
            --assembly {input.a} \
            --bin_motifs {input.b} \
            --contig_bins {input.c} \
            --write_bins \
            --out {OUTDIR}/{wildcards.sample}/nanomotif_binnary/
        """

rule checkm2_decontaminate:
    conda: "checkm2"
    input:
        os.path.join(OUTDIR,"{sample}","nanomotif_binnary",  "bin_contamination.tsv"),
    output:
        os.path.join(OUTDIR, "{sample}", "nanomotif_binnary", "checkm2_decon", "quality_report.tsv")
    resources:
        mem = "40G",
        walltime = "08:00:00",
        nodetype = config["partition"],
    threads: 20
    shell:
        """
            export CHECKM2DB="data/databases/CheckM2_database/uniref100.KO.1.dmnd"
            checkm2 predict \
                --threads {threads} \
                --input {OUTDIR}/{wildcards.sample}/nanomotif_binnary/detect_contamination_bins/ \
                --force \
                -x ".fa" \
                --output-directory {OUTDIR}/{wildcards.sample}/nanomotif_binnary/checkm2_decon/
        """
        

rule include:
    conda:
        "nanomotif_dev"
    input:
        a = lambda wildcards: config["samples"][wildcards.sample]["assembly"],
        p = lambda wildcards: config["samples"][wildcards.sample]["pileup"],
        b = lambda wildcards: config["samples"][wildcards.sample]["bin_motifs"],
        bin_con = os.path.join(OUTDIR,"{sample}","nanomotif_binnary",  "bin_contamination.tsv"),
        c = lambda wildcards: config["samples"][wildcards.sample]["contig_bin"],
        m = os.path.join(OUTDIR,"{sample}","nanomotif_binnary", "motifs-read-methylation-copied.done")
    output:
        os.path.join(OUTDIR,"{sample}","nanomotif_binnary",  "include_contigs.tsv"),
    threads:
        5
    resources:
        mem = "40G",
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
            --out {OUTDIR}/{wildcards.sample}/nanomotif_binnary/
        """



