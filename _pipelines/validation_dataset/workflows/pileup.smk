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
            os.path.join(OUTDIR, "{dataset}","1_cov.bam_0_data_cov.txt"),
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

rule coverage_contigs_bed:
    conda:
        "envs/mapreads.yaml"
    input: 
        bam = os.path.join(OUTDIR, "{dataset}", "1_cov.bam"),
    output:
        temp(os.path.join(OUTDIR, "{dataset}","contigs.bed"))
    threads:
        1
    resources:
        mem = "1G",
        nodetype = "general",
        walltime = "10:00",
    shell:
        """
        samtools view -H {input} | grep "@SQ" | awk -F'\\t' '{{split($2,a,":"); split($3,b,":"); print a[2]"\\t0\\t"b[2]}}' > {output}
        """
        
rule coverage_extract:
    conda:
        "envs/bedtools.yaml",
    input:
        contigs = os.path.join(OUTDIR, "{dataset}", "contigs.bed"),
        bam = os.path.join(OUTDIR, "{dataset}", "1_cov.bam"),
    output:
        os.path.join(OUTDIR, "{dataset}","1_cov.bam_0_data_cov.txt"),
    threads: 1
    resources:
        mem = "60G",
        nodetype = "high-mem",
        walltime = "9:00:00",
    shell:
        """
        bedtools coverage -a {input.contigs} -b {input.bam} -mean > {output}
        """

