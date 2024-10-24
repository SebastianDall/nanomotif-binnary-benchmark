import os
import glob
import datetime

SNAKEDIR = os.path.dirname(workflow.snakefile)
configfile: os.path.join(SNAKEDIR, "..", "config", "config-subsample.yaml")


OUTDIR = config["outdir"]
VAL_OUTDIR = config["validation_outdir"]
SAMPLES = config["samples"]


rule all:
    input:
        expand(os.path.join(VAL_OUTDIR, "nanoq", "{sample}.json"), sample=SAMPLES.keys()),
        expand(os.path.join(VAL_OUTDIR, "qc", "{sample}_qc_min_l_" + str(config["min_l"]) + "_min_q_" + str(config["min_q"]) + ".fastq"), sample=SAMPLES.keys()),
        expand(os.path.join(VAL_OUTDIR, "nanoq", "{sample}_qc_min_l_" + str(config["min_l"]) + "_min_q_" + str(config["min_q"]) + ".json"), sample=SAMPLES.keys())

def get_fastq(wildcards):
    sample = SAMPLES[wildcards.sample]["fastq"]
    
    return sample

rule nanoq:
    input:
        i = get_fastq
    output:
        o = os.path.join(VAL_OUTDIR, "nanoq", "{sample}.json")
    threads: 1
    resources:
        mem = "15G",
        walltime = "1:00:00",
        nodetype = "general"
    params:
        nanoq = os.path.join(SNAKEDIR, "envs", "nanoq")
    shell:
        """
        {params.nanoq} -i {input.i} --json -f -s > {output.o}
        """


rule trim:
    input:
        i = get_fastq
    output:
        o = os.path.join(VAL_OUTDIR, "qc", "{sample}_qc_min_l_" + str(config["min_l"]) + "_min_q_" + str(config["min_q"]) + ".fastq"),
        L = os.path.join(VAL_OUTDIR, "qc", "{sample}_qc_min_l_" + str(config["min_l"]) + "_min_q_" + str(config["min_q"]) + ".txt")
    threads: 1
    resources:
        mem = "15G",
        walltime = "1:00:00",
        nodetype = "general"
    params:
        nanoq = os.path.join(SNAKEDIR, "envs", "nanoq"),
        min_length = config["min_l"],
        min_q = config["min_q"],
    shell:
        """
        {params.nanoq} -i {input.i} -l {params.min_length} -q {params.min_q} -L {output.L} > {output.o}
        """


rule nanoq_trim:
    input:
        i = os.path.join(VAL_OUTDIR, "qc", "{sample}_qc_min_l_" + str(config["min_l"]) + "_min_q_" + str(config["min_q"]) + ".fastq")
    output:
        o = os.path.join(VAL_OUTDIR, "nanoq", "{sample}_qc_min_l_" + str(config["min_l"]) + "_min_q_" + str(config["min_q"]) + ".json")
    threads: 1
    resources:
        mem = "15G",
        walltime = "1:00:00",
        nodetype = "general"
    params:
        nanoq = os.path.join(SNAKEDIR, "envs", "nanoq")
    shell:
        """
        {params.nanoq} -i {input.i} --json -f -s > {output.o}
        """