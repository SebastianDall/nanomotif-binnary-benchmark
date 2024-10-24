import os
import glob


SNAKEDIR = os.path.dirname(workflow.snakefile)
configfile: os.path.join(SNAKEDIR, "..", "config", "config-subsample.yaml")

OUTDIR = config["outdir"]
SAMPLES = config["samples"]

DATASETS = config["dataset"]


def input():
    samples = []
    for dataset in DATASETS:
        for sample in DATASETS[dataset]["samples"]:
            samples.append(os.path.join(OUTDIR, dataset, "contig_mapping", sample, "coverage.csv"))
    return samples


rule all:
    input:
        input(),
        expand(os.path.join(OUTDIR, "{dataset}", "contig_mapping", "mapped_contig_bin.tsv"), dataset=DATASETS.keys())


rule minimap2:
    conda:
        "envs/mapreads.yaml"
    input: 
        ref = lambda wildcards: DATASETS[wildcards.dataset]['ref'],
        query = lambda wildcards: SAMPLES[wildcards.sample]['assembly']
    output:
        os.path.join(OUTDIR, "{dataset}", "contig_mapping", "{sample}","mapping.sam")
    threads: 10,
    resources:
        mem = "30G",
        walltime = "7:00:00",
        nodetype = "general"
    shell:
        """
        minimap2 -a {input.ref} {input.query} > {output}
        """


rule samtools_bam:
    conda:
        "envs/mapreads.yaml"
    input:
        sam = os.path.join(OUTDIR, "{dataset}", "contig_mapping", "{sample}","mapping.sam")
    output:
        os.path.join(OUTDIR, "{dataset}", "contig_mapping", "{sample}", "mapping.bam")
    threads:
        5
    resources:
        mem = "10G",
        walltime = "7:00:00",
        nodetype = "general"
    shell:
        """
        samtools view -bS {input.sam} > {output}
        """

rule samtools_sort:
    conda:
        "envs/mapreads.yaml"
    input:
        bam = os.path.join(OUTDIR, "{dataset}", "contig_mapping", "{sample}", "mapping.bam")
    output:
        os.path.join(OUTDIR, "{dataset}", "contig_mapping", "{sample}", "mapping.sorted.bam")
    threads: 5
    resources:
        mem = "10G",
        walltime = "7:00:00",
        nodetype = "general"
    shell:
        """
        samtools sort {input.bam} -o {output}
        """

rule samtools_index:
    conda:
        "envs/mapreads.yaml"
    input:
        bam = os.path.join(OUTDIR, "{dataset}", "contig_mapping", "{sample}","mapping.sorted.bam")
    output:
        os.path.join(OUTDIR, "{dataset}", "contig_mapping", "{sample}","mapping.sorted.bam.bai")
    threads: 5
    resources:
        mem = "10G",
        walltime = "7:00:00",
        nodetype = "general"
    shell:
        """
        samtools index {input.bam}
        """
    
rule samtools_hq_map:
    conda:
        "envs/mapreads.yaml"
    input:
        bam = os.path.join(OUTDIR, "{dataset}", "contig_mapping", "{sample}","mapping.sorted.bam"),
        bai = os.path.join(OUTDIR, "{dataset}", "contig_mapping", "{sample}","mapping.sorted.bam.bai")
    output:
        os.path.join(OUTDIR, "{dataset}", "contig_mapping", "{sample}","mapping.hq.bam")
    threads: 5
    resources:
        mem = "10G",
        walltime = "7:00:00",
        nodetype = "general"
    shell:
        """
        samtools view -b -q 30 {input.bam} > {output}
        """



rule bedtools_coverage:
    conda:
        "envs/bedtools.yaml"
    input:
        bam = os.path.join(OUTDIR, "{dataset}", "contig_mapping", "{sample}","mapping.hq.bam"),
    output:
        os.path.join(OUTDIR, "{dataset}", "contig_mapping", "{sample}","coverage.bed")
    threads: 5
    resources:
        mem = "10G",
        walltime = "7:00:00",
        nodetype = "general"
    shell:
        """
        bedtools genomecov -ibam {input.bam} -bg > {output}
        """

rule filter_coverage:
    input:
        bed = os.path.join(OUTDIR, "{dataset}", "contig_mapping", "{sample}","coverage.bed"),
        info = os.path.join(OUTDIR, "{dataset}", "flye", "assembly_info.txt")
    output:
        csv = os.path.join(OUTDIR, "{dataset}", "contig_mapping", "{sample}","coverage.csv")
    params:
        coverage_mapper = config["coverage_mapper"],
        coverage_threshold = config["coverage_threshold"]
    threads:
        5
    resources:
        mem = "10G",
        walltime = "7:00:00",
        nodetype = "general"
    shell:
        """
        {params.coverage_mapper} -l {input.info} -c {input.bed} --alignment-pct {params.coverage_threshold} --name {wildcards.sample} --out {output}
        """

def getCoverageSamples(wildcards):
    return glob.glob(os.path.join(OUTDIR, wildcards.dataset, "contig_mapping", "*", "coverage.csv"))

rule gather_coverage:
    input:
        getCoverageSamples
    output:
        os.path.join(OUTDIR, "{dataset}", "contig_mapping", "coverage.csv")
    shell:
        """
        echo "contig,coverage,sample" > {output}
        cat {input} >> {output}
        """


rule mapped_contig_bin:
    input:
        os.path.join(OUTDIR, "{dataset}", "contig_mapping", "coverage.csv")        
    output:
        os.path.join(OUTDIR, "{dataset}", "contig_mapping", "mapped_contig_bin.tsv")
    threads:
        1
    run:
        import pandas as pd
        contig_mapping = pd.read_csv(input[0], sep = ",")
        # Filter rows where coverage > 0.9
        filtered_data = contig_mapping[contig_mapping['coverage'] > 0.9]

        # Group by 'contig' and keep groups with only one entry
        grouped_data = filtered_data.groupby('contig').filter(lambda x: len(x) == 1)

        # Rename 'sample' column to 'bin'
        grouped_data = grouped_data.rename(columns={'sample': 'bin'})

        # Drop the 'coverage' column
        grouped_data = grouped_data.drop(columns=['coverage'])

        # Sort by 'bin'
        sorted_data = grouped_data.sort_values(by='bin')

        # Write the DataFrame to a TSV file
        sorted_data.to_csv(output[0], sep='\t', index=False)
