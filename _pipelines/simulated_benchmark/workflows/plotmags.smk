import os
import glob
import datetime

SNAKEDIR = os.path.dirname(workflow.snakefile)
configfile: os.path.join(SNAKEDIR, "..", "config", "config.yaml")

BENCHMARK_DIR = config["output_dir"]
BENCHMARKS = [config["benchmarks"]]
sample_path = os.path.join(BENCHMARK_DIR, BENCHMARKS[0], "*")

SAMPLES = glob.glob(sample_path)
SAMPLES = [os.path.basename(s) for s in SAMPLES]


rule all:
    input:
        expand(os.path.join(BENCHMARK_DIR, "{benchmark}", "{sample}", "MAGs", "done.txt"), sample=SAMPLES, benchmark=BENCHMARKS),


def find_genomad(wildcards):
    files = glob.glob(f'{os.path.join(config["data_dir"], wildcards.sample, "genomad")}/**/*_aggregated_classification.tsv', recursive=True)
    # remove files that contain provirus
    files = [f for f in files if "provirus" not in f]
    return files

rule plot_MAGs:
    input:
        motifs_scored = lambda wildcards: os.path.join(config["general_data_dir"], "Nanomotif", wildcards.sample, "motifs-scored.tsv"),
        bin_motifs = lambda wildcards: os.path.join(config["general_data_dir"], "Nanomotif", wildcards.sample, "bin-motifs.tsv"),
        checkm2 = lambda wildcards: os.path.join(BENCHMARK_DIR, wildcards.benchmark, wildcards.sample, "checkm2", "quality_report.tsv"),
        assembly = lambda wildcards: config['samples'][wildcards.sample]['assembly'],
        cov = lambda wildcards: os.path.join(config["general_data_dir"], "SemiBin2", wildcards.sample, "1_cov.bam_0_data_cov.csv"),
        contig_bins = lambda wildcards: os.path.join(BENCHMARK_DIR, wildcards.benchmark, wildcards.sample, "contig_bins.tsv"),
        genomad = find_genomad,
    output:
        touch(os.path.join(BENCHMARK_DIR, "{benchmark}", "{sample}", "MAGs", "done.txt"))
    threads:
        3
    resources:
        mem = "15G",
        walltime = "00:25:00",
        nodetype = config["partition_general"],
    params:
        script = os.path.join(SNAKEDIR, "src","plot_MAGs.R"),
    shell:
        """
        Rscript {params.script} --motifs_scored {input.motifs_scored} \
            --bin_motifs {input.bin_motifs} \
            --checkm2 {input.checkm2} \
            --contig_bins {input.contig_bins} \
            --assembly {input.assembly} \
            --cov {input.cov} \
            --genomad {input.genomad} \
            --output {BENCHMARK_DIR}/{wildcards.benchmark}/{wildcards.sample}/MAGs  
        """