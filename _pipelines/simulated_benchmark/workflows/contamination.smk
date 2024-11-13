import os
import glob
import datetime

SNAKEDIR = os.path.dirname(workflow.snakefile)
configfile: os.path.join(SNAKEDIR, "..", "config", "config-contamination.yaml")

# Extract sample names from the YAML configuration
# BENCHMARKS = [benchmark for benchmark in config['benchmarks'].keys()]
DATETIME = datetime.datetime.now().strftime('%Y%m%d_%H%M')
BASELINE = config["baseline_dir"]
OUTDIR = config["outdir"]


def benchmark_inputs():
    inputs = []

    for sample in config["samples"].keys():
        for benchmark in config["samples"][sample]["benchmarks"].keys():
            inputs.append(os.path.join(BASELINE, "detect_contamination", sample, benchmark, "bin_contamination.tsv"))
            inputs.append(os.path.join(OUTDIR, sample, benchmark, "bin_contamination.tsv"))
    return inputs


rule all:
    input:
        benchmark_inputs(),
        expand(os.path.join(BASELINE, "motif_discovery", "{sample}",  "original_contig_bin","motifs-scored.tsv"),
               sample = config["samples"].keys()),
        # expand(os.path.join(BASELINE, "detect_contamination", "{sample}", "{benchmark}", "bin_contamination.tsv"),
        #        sample = config["samples"].keys(),
        #        benchmark = config["samples"][sample]["benchmarks"])
        # expand(os.path.join(BASELINE, "detect_contamination", "{benchmark}", "bin_contamination.tsv"),
        #        benchmark = BENCHMARKS),
        # expand(os.path.join(OUTDIR, "{benchmark}", "bin_contamination.tsv"),
        #        benchmark = BENCHMARKS)

       

rule nanomotif_motif_discovery:
    conda:
        "envs/nanomotif.yaml"
    input:
        a = lambda wildcards: config["samples"][wildcards.sample]["assembly"],
        p = lambda wildcards: config["samples"][wildcards.sample]["pileup"],
        bins = lambda wildcards: config["samples"][wildcards.sample]["benchmarks"]["original_contig_bin"],
    output:
        o = os.path.join(BASELINE, "motif_discovery", "{sample}", "original_contig_bin", "bin-motifs.tsv"),
        motifs = os.path.join(BASELINE, "motif_discovery", "{sample}", "original_contig_bin","motifs.tsv"),
        m = os.path.join(BASELINE, "motif_discovery", "{sample}", "original_contig_bin","motifs-scored.tsv"),
    threads:
        40
    resources:
        mem = "40G",
        walltime = "1-00:00:00",
        nodetype = config["partition"],
    shell:
        """
        nanomotif motif_discovery \
            -t {threads} \
            --out {BASELINE}/motif_discovery/{wildcards.sample}/original_contig_bin/ \
            {input.a} {input.p} {input.bins}
        """

rule copy_contig_bin_baseline:
    input:
        lambda wildcards: config["samples"][wildcards.sample]["benchmarks"][wildcards.benchmark],
    output:
        os.path.join(BASELINE, "detect_contamination", "{sample}", "{benchmark}", "contig_bin.tsv"),
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
        


rule nanomotif_contamination:
    conda:
        "envs/nanomotif.yaml"
    input:
        m = os.path.join(BASELINE, "motif_discovery", "{sample}", "original_contig_bin", "motifs-scored.tsv"),
        b = os.path.join(BASELINE, "motif_discovery", "{sample}", "original_contig_bin", "bin-motifs.tsv"),
        c = os.path.join(BASELINE, "detect_contamination", "{sample}", "{benchmark}", "contig_bin.tsv"),
    output:
        os.path.join(BASELINE, "detect_contamination", "{sample}", "{benchmark}", "bin_contamination.tsv"),
    threads:
        20
    resources:
        mem = "40G",
        walltime = "08:00:00",
        nodetype = config["partition"],
    params:
    shell:
        """
        nanomotif detect_contamination \
            --motifs_scored {input.m} \
            --bin_motifs {input.b} \
            --contig_bins {input.c} \
            --out {BASELINE}/detect_contamination/{wildcards.sample}/{wildcards.benchmark}         
        """

rule nanomotif_contamination_dev:
    conda:
        "nanomotif_dev"
    input:
        a = lambda wildcards: config["samples"][wildcards.sample]["assembly"],
        p = lambda wildcards: config["samples"][wildcards.sample]["pileup"],
        b = os.path.join(BASELINE, "motif_discovery", "{sample}", "original_contig_bin", "bin-motifs.tsv"),
        c = os.path.join(OUTDIR,"{sample}","{benchmark}",  "contig_bin.tsv"),
    output:
        os.path.join(OUTDIR,"{sample}","{benchmark}",  "bin_contamination.tsv"),
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
        nanomotif detect_contamination \
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
