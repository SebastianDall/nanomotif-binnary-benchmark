import os
import glob
import datetime
import pandas as pd

SNAKEDIR = os.path.dirname(workflow.snakefile)
configfile: os.path.join(SNAKEDIR, "..", "config", "config-subsample.yaml")

OUTDIR = config["outdir"]
VAL_OUTDIR = config["validation_outdir"]
SAMPLES = config["samples"]

DATASETS = config["dataset"]

def define_subsamples():
    subsamples = []
    for dataset in DATASETS:
        for sample in DATASETS[dataset]["samples"]:
            p = os.path.join(VAL_OUTDIR, dataset, sample + "_subsampled.fastq")
            subsamples.append(p)
    return subsamples
print(list(DATASETS.keys()))
rule all:
    input:
        define_subsamples(),
        expand(os.path.join(OUTDIR, "{dataset}", "concatenated.fastq.gz"), dataset=DATASETS.keys()),
        expand(os.path.join(OUTDIR, "{dataset}", "assembly.fasta"), dataset=DATASETS.keys()),


num_reads_dict = {}
for dataset in DATASETS:
    rank_abundance = pd.read_csv(DATASETS[dataset]["rank_abundance_file"])
    num_reads_dict[dataset] = {}

    for sample in DATASETS[dataset]["samples"]:
        # save num reads for each sample for each dataset
        num_reads_required = int(rank_abundance[rank_abundance["id"] == sample]["num_reads_required"].values[0])
        num_reads = int(rank_abundance[rank_abundance["id"] == sample]["num_reads"].values[0])
        proportion = min((num_reads_required / num_reads) + 0.1, 1)

        # Initialize the sample dictionary
        num_reads_dict[dataset][sample] = {}

        num_reads_dict[dataset][sample]["num_reads"] = num_reads_required
        num_reads_dict[dataset][sample]["p"] = proportion

rule subsample_fastq:
    conda: "envs/seqkit.yaml"
    input:
        fq = lambda wildcards: SAMPLES[wildcards.sample]["fastq"]
    output:
        os.path.join(VAL_OUTDIR,  "{dataset}", "{sample}_subsampled.fastq")
    threads: 10
    resources:
        mem = "30G",
        walltime = "7:00:00",
        nodetype = "general"
    params:
        num_reads = lambda wildcards: num_reads_dict[wildcards.dataset][wildcards.sample]["num_reads"],
        p = lambda wildcards: num_reads_dict[wildcards.dataset][wildcards.sample]["p"]
    shell:
        """
        # Subsampling using seqkit
        seqkit sample -p {params.p} {input.fq} -s 11 | seqkit sample -n {params.num_reads} -s 11 > {output}
        """


rule check_num_reads:
    input:
        os.path.join(VAL_OUTDIR,  "{dataset}", "{sample}_subsampled.fastq")
    output:
        touch(os.path.join(VAL_OUTDIR,  "{dataset}", "{sample}_subsampled.fastq.num_reads_ok"))
    threads: 1
    params:
        expected_num_reads = lambda wildcards: num_reads_dict[wildcards.dataset][wildcards.sample]["num_reads"]
    resources:
        mem = "1G",
        walltime = "1:00:00",
        nodetype = "general"
    run:
        num_reads = 0
        with open(input[0]) as f:
            for line in f:
                if line.startswith("@"):
                    num_reads += 1
        if num_reads >= params.expected_num_reads - 10:
            shell(f"echo {num_reads} > {output[0]}")
            
        else:
            raise ValueError(f"Number of reads in {input[0]} is {num_reads}, expected {params.expected_num_reads}")


rule concatenate_fastq:
    input:
        fastqs = lambda wildcards: expand(
            os.path.join(VAL_OUTDIR,  "{dataset}", "{sample}_subsampled.fastq"), 
            dataset=wildcards.dataset, sample=DATASETS[wildcards.dataset]["samples"]
        ),
        p = lambda wildcards: expand(
            os.path.join(VAL_OUTDIR,  "{dataset}", "{sample}_subsampled.fastq.num_reads_ok"),
            dataset=wildcards.dataset, sample=DATASETS[wildcards.dataset]["samples"]
        )
    output:
        os.path.join(OUTDIR, "{dataset}", "concatenated.fastq")
    resources:
        mem = "15G",
        walltime = "7:00:00",
        nodetype = "general"
    shell:
        """
        # Concatenate and gzip the fastq files
        cat {input.fastqs} > {output}
        """



rule assembly:
    conda:
        "envs/flye.yaml"
    input:
        os.path.join(OUTDIR, "{dataset}", "concatenated.fastq")
    output:
        # directory = 
        a = os.path.join(OUTDIR, "{dataset}", "assembly.fasta")
    threads: 40
    params:
        directory = lambda wildcards: os.path.join(OUTDIR, wildcards.dataset)
    resources:
        mem = "100G",
        walltime = "2-00:00:00",
        nodetype = "general"
    shell:
        """
        flye --nano-hq {input} --out-dir {params.directory} --threads {threads} \
        --meta
        """



rule minimap2:
    conda:
        "envs/mapreads.yaml"
    input: 
        ref = os.path.join(OUTDIR, "{dataset}", "assembly.fasta"),
        query = os.path.join(VAL_OUTDIR, "renamed_assemblies", "{sample}.fasta")
    output:
        "data/alignments/{dataset}/{sample}.sam"
    threads: 20
    resources:
        mem = "40",
        walltime = "2-00:00:00",
        nodetype = "general"
    shell:
        """
        minimap2 -a {input.ref} {input.query} > {output}
        """