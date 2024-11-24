import os
import glob
import datetime
import polars as pl

SNAKEDIR = os.path.dirname(workflow.snakefile)
configfile: os.path.join(SNAKEDIR, "..", "config", "config-contamination.yaml")


# Extract sample names from the YAML configuration
# BENCHMARKS = [benchmark for benchmark in config['benchmarks'].keys()]
DATETIME = datetime.datetime.now().strftime('%Y%m%d_%H%M')
BASELINE = config["baseline_dir"]
OUTDIR = config["outdir"]

# include: "rules/motif_discovery.smk"

def benchmark_inputs():
    inputs = []

    for sample in config["samples"].keys():
        for benchmark in config["samples"][sample]["benchmarks"].keys():
            inputs.append(os.path.join(OUTDIR, sample, benchmark, "bin_contamination.tsv"))
    return inputs


rule all:
    input:
        expand(
            os.path.join(OUTDIR, "{sample}", "benchmarks.done"),
            sample=config["samples"].keys()
        )
        
        # benchmark_inputs(),
       
        


rule create_contamination_pileup:
    input:
        p = lambda wildcards: config["samples"][wildcards.sample]["pileup"]
    output:
        p = os.path.join(BASELINE, "contamination_files", "{sample}", "contamination_pileup.bed")
    threads:
        10
    resources:
        mem = "140G",
        walltime = "01:00:00",
        nodetype = "high-mem",
    run:
        import polars as pl

        pileup = pl.read_csv(input[0], separator="\t", has_header = False)
        con_pileup = pileup.with_columns(
            (pl.lit("contamination_") + pl.col("column_1")).alias("column_1")
        )

        p_con = pl.concat([pileup, con_pileup])
        p_con.write_csv(output[0], separator="\t", include_header = False)



rule create_contamination_assembly:
    input:
        a = lambda wildcards: config["samples"][wildcards.sample]["assembly"]
    output:
        a = os.path.join(BASELINE, "contamination_files", "{sample}", "contamination_assembly.fasta")
    threads:
        10
    resources:
        mem = "15G",
        walltime = "01:00:00",
        nodetype = config["partition"],
    run:
        from Bio import SeqIO
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq

        def create_contamination_assembly(input_file, output_file):
            original_contigs = list(SeqIO.parse(input_file, "fasta"))

            duplicated_contigs = []
            for contig in original_contigs:
                duplicated_contig = SeqRecord(
                    Seq(str(contig.seq)),
                    id = f"contamination_{contig.id}",
                    description = contig.description
                )
                duplicated_contigs.append(duplicated_contig)

            combined_contigs = original_contigs + duplicated_contigs

            SeqIO.write(combined_contigs, output_file, "fasta")

        create_contamination_assembly(input[0], output[0])

rule create_motifs_scored_file:
    input:
        a = os.path.join(BASELINE, "contamination_files", "{sample}", "contamination_assembly.fasta"),
        p = os.path.join(BASELINE, "contamination_files", "{sample}", "contamination_pileup.bed"),
        b = lambda wildcards: config["samples"][wildcards.sample]["bin_motifs"],
    output:
        m = os.path.join(BASELINE, "contamination_files", "{sample}", f"motifs-scored-read-methylation-{config['min_valid_read_cov']}.tsv")
    threads:
        40
    resources:
        mem = "200G",
        walltime = "08:00:00",
        nodetype = "high-mem",
    params:
        MU = os.path.join(SNAKEDIR, "src", "methylation_utils")
    run:
        from methylation_utils_wrapper.utils import run_methylation_utils
        import sys

        motifs = pl.read_csv(input[2], separator = "\t")
        motifs = motifs\
            .with_columns((pl.col("motif") + "_" + pl.col("mod_type") + "_" + pl.col("mod_position").cast(pl.Utf8)).alias("motif_mod"))\
            .with_columns([
                (pl.col("n_mod_bin") + pl.col("n_nomod_bin")).alias("n_motifs"),
                (pl.col("n_mod_bin") / (pl.col("n_mod_bin") + pl.col("n_nomod_bin"))).alias("mean_methylation")
            ])\
            .filter(pl.col("mean_methylation") >= config["mean_methylation_cutoff"])\
            .filter(pl.col("n_motifs") >= config["n_motif_bin_cutoff"])\
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
        m = os.path.join(BASELINE, "contamination_files", "{sample}", f"motifs-scored-read-methylation-{config['min_valid_read_cov']}.tsv")
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
        m = os.path.join(BASELINE, "contamination_files", "{sample}", f"motifs-scored-read-methylation-{config['min_valid_read_cov']}.tsv"),
        c = lambda wildcards: config["samples"][wildcards.sample]["contig_bin"]
    output:
        directory(os.path.join(BASELINE, "benchmark_files", "{sample}", "contamination"))
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
            python {params.PY} {input.c} {input.m} {params.n_benchmarks} contamination {BASELINE}/"benchmark_files"/{wildcards.sample}/
        """

        

    
        
rule symlink_contig_bin_dev:
    input:
        os.path.join(BASELINE, "benchmark_files", "{sample}", "contamination", "{benchmark}.tsv")
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
        

rule nanomotif_contamination_dev:
    conda:
        "nanomotif_dev"
    input:
        a = os.path.join(BASELINE, "contamination_files", "{sample}", "contamination_assembly.fasta"),
        p = os.path.join(BASELINE, "contamination_files", "{sample}", "contamination_pileup.bed"),
        b = os.path.join(BASELINE, "motif_discovery", "{sample}", "original_contig_bin", "bin-motifs.tsv"),
        c = os.path.join(OUTDIR,"{sample}","{benchmark}",  "contig_bin.tsv"),
        m = os.path.join(OUTDIR,"{sample}","{benchmark}", "motifs-read-methylation-copied.done")
    output:
        os.path.join(OUTDIR,"{sample}","{benchmark}",  "bin_contamination.tsv"),
    threads:
        10
    resources:
        mem = "20G",
        walltime = "08:00:00",
        nodetype = "default-op",
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

def find_benchmarks(wildcards):
    cp = checkpoints.create_benchmark_contig_bin_files.get(**wildcards).output[0]

    return expand(
        os.path.join(OUTDIR, "{sample}", "{benchmark}", "bin_contamination.tsv"),
        sample = wildcards.sample, benchmark = glob_wildcards(os.path.join(cp, "{benchmark}.tsv")).benchmark
    )

rule aggregate_benchmarks:
    input:
        find_benchmarks
    output:
        touch(os.path.join(OUTDIR, "{sample}", "benchmarks.done"))
    threads:
        1
    resources:
        mem = "2G",
        walltime = "00:02:00",
        nodetype = "default-op",
