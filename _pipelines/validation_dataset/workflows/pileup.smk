import os
import glob
import datetime

SNAKEDIR = os.path.dirname(workflow.snakefile)
configfile: os.path.join(SNAKEDIR, "..", "config", "config.yaml")


VALIDATION_DATASETS = config["validation_datasets"]
OUTDIR = config["outdir"]
VAL_OUTDIR = config["validation_outdir"]

rule all:
    input:
        expand(
            os.path.join(OUTDIR, "{dataset}", "mod_pileup.bed"),
            dataset=VALIDATION_DATASETS.keys()
        ),
        expand(
            os.path.join(OUTDIR, "{dataset}", "1_cov.bam"),
            dataset=VALIDATION_DATASETS.keys()
        )


rule map_mod_calls:
    conda: 
        "envs/mapreads.yaml"
    input:
        assembly = lambda wildcards: VALIDATION_DATASETS[wildcards.dataset]["assembly"],
        reads = lambda wildcards: VALIDATION_DATASETS[wildcards.dataset]["fastq"]
    output:
        bam = os.path.join(OUTDIR, "{dataset}", "mod.bam")
    threads: 40,
    params:
        mem = "50g",
        ref = "10g"
    resources:
        mem ="400G",
        nodetype="high-mem",
        walltime="2-00:00:00",
    shell:
        """
        minimap2 -ax map-ont -t {threads} -I {params.mem} -K {params.ref} -y {input.assembly} {input.reads} | \
        samtools view -bS | \
        samtools sort > {output.bam}
        """


rule index_bam:
    conda: 
        "envs/mapreads.yaml"
    input:
        bam = os.path.join(OUTDIR, "{dataset}", "mod.bam")
    output:
        bai = os.path.join(OUTDIR, "{dataset}", "mod.bam.bai")
    threads: 10,
    resources:
        mem ="10G",
        nodetype="general",
        walltime="1-00:00:00",
    shell:
        """
        samtools index -@ {threads} {input.bam}
        """

rule pileup:
    conda:
        "modkit"
    input:
        bam = os.path.join(OUTDIR, "{dataset}", "mod.bam"),
        bai = os.path.join(OUTDIR, "{dataset}", "mod.bam.bai"),
    output:
        pileup = os.path.join(OUTDIR, "{dataset}", "mod_pileup.bed"),
    threads: 40,
    resources:
        mem ="400G",
        nodetype="high-mem",
        walltime="3-00:00:00",
    shell:
        """
        modkit pileup -t {threads} --only-tabs {input.bam} {output.pileup} 
        """

rule map_cov:
    conda: 
        "envs/mapreads.yaml"
    input:
        assembly = lambda wildcards: VALIDATION_DATASETS[wildcards.dataset]["assembly"],
        reads = lambda wildcards: VALIDATION_DATASETS[wildcards.dataset]["fastq"]
    output:
        bam = os.path.join(OUTDIR, "{dataset}", "1_cov.bam"),
    threads: 40,
    params:
        mem = "50g",
        ref = "10g"
    resources:
        mem ="400G",
        nodetype="high-mem",
        walltime="2-00:00:00",
    shell:
        """
        minimap2 -ax map-ont -t {threads} -I {params.mem} -K {params.ref} -y {input.assembly} {input.reads} | \
        samtools view -bS -F 2308 | \
        samtools sort > {output.bam}
        """