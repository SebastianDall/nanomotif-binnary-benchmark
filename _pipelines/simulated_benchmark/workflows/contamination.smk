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

include: "rules/motif_discovery.smk"

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


rule create_contamination_pileup:
    input:
        p = lambda wildcards: config["samples"][wildcards.sample]["pileup"]
    output:
        p = os.path.join(BASELINE, "contamination_files", "{sample}", "contamination_pileup.bed")
    threads:
        10
    resources:
        mem = "65G",
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

rule nanomotif_contamination_dev:
    conda:
        "nanomotif_dev"
    input:
        a = os.path.join(BASELINE, "contamination_files", "{sample}", "contamination_assembly.fasta"),
        p = os.path.join(BASELINE, "contamination_files", "{sample}", "contamination_pileup.bed"),
        b = os.path.join(BASELINE, "motif_discovery", "{sample}", "original_contig_bin", "bin-motifs.tsv"),
        c = os.path.join(OUTDIR,"{sample}","{benchmark}",  "contig_bin.tsv"),
    output:
        os.path.join(OUTDIR,"{sample}","{benchmark}",  "bin_contamination.tsv"),
    threads:
        40
    resources:
        mem = "100G",
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
