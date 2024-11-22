import os
import glob
import datetime
import pandas as pd

SNAKEDIR = os.path.dirname(workflow.snakefile)
configfile: os.path.join(SNAKEDIR, "..", "config", "config-fragmentation_benchmark.yaml")

OUTDIR = config["outdir"]


rule all:
  input:
    os.path.join(OUTDIR, "concatenated_assembly.fasta"),
    os.path.join(OUTDIR, "concatenated_reads.fastq"),
    expand(os.path.join(OUTDIR, "contig_bin_{fragment_size}kb.tsv"),fragment_size = config["fragment_size"]),
    expand(
      os.path.join(OUTDIR, "subset_reads/{coverage}x/{sample}_{coverage}x.fastq"),
      coverage = config["coverage"],
      sample = config["samples"]
      
    ),
    expand(
      os.path.join(OUTDIR,"chunked_assembly/{fragment_size}kb/{sample}-{fragment_size}kb.fasta"),
      fragment_size = config["fragment_size"],
      sample = config["samples"]
      
    )

################################################################################
# Monoculture data subsetting
 
rule chunk_reference_to_fragment_size:
  conda: "envs/seqkit.yaml",
  input: lambda wildcards: config["samples"][wildcards.sample]["assembly"]
  output: os.path.join(OUTDIR,"chunked_assembly/{fragment_size}kb/{sample}-{fragment_size}kb.fasta")
  threads: 1
  resources:
      mem = "5G",
      walltime = "00:30:00",
      nodetype = config["nodetype"]
  shell:
      """
      seqkit sliding -s {wildcards.fragment_size}000 -W {wildcards.fragment_size}000 {input} \
      | sed '/^>/ s/:/_/g' \
      | sed '/^>/ s/^>/>{wildcards.sample}_/' > {output}
      """

rule subset_reads_to_x_coverage:
  conda: "envs/mapreads.yaml",
  input:
      reads=lambda wildcards: config["samples"][wildcards.sample]["bam"],
      coverage=lambda wildcards: config["samples"][wildcards.sample]["coverage"]
  output: os.path.join(OUTDIR,"subset_reads/{coverage}x/{sample}_{coverage}x.bam")
  threads: 1
  resources:
      mem = "5G",
      walltime = "00:30:00",
      nodetype = config["nodetype"]
  shell:
      """
      genome_coverage=$(awk 'NR==3 {{print $3}}' {input.coverage})
      subset_fraction=$(echo "scale=2; {wildcards.coverage} / $genome_coverage" | bc)
      if [ $subset_fraction -gt 1 ]; then subset_fraction=1; fi
      samtools view -s $subset_fraction -b {input.reads} > {output}
      """
 
rule subset_reads_bam_to_fastq:
  conda: "envs/mapreads.yaml",
  input: os.path.join(OUTDIR,"subset_reads/{coverage}x/{sample}_{coverage}x.bam")
  output: os.path.join(OUTDIR,"subset_reads/{coverage}x/{sample}_{coverage}x.fastq")
  threads: 1
  resources:
      mem = "5G",
      walltime = "00:30:00",
      nodetype = config["nodetype"]
  shell:
      """
      samtools fastq -T MM,ML {input} > {output}
      """
 
rule concatenate_reads:
  input:
    expand(
      os.path.join(OUTDIR, "subset_reads/{coverage}x/{sample}_{coverage}x.fastq"),
      coverage = config["coverage"],
      sample = config["samples"].keys()
    ),
  output:
    os.path.join(OUTDIR, "concatenated_reads.fastq")
  threads: 1
  resources:
      mem = "5G",
      walltime = "10:30:00",
      nodetype = config["nodetype"]
  shell:
    """
      cat {input} > {output}
    """

rule concatenate_assembly:
  input:
    expand(
      os.path.join(OUTDIR,"chunked_assembly/{fragment_size}kb/{sample}-{fragment_size}kb.fasta"),
      fragment_size = config["fragment_size"],
      sample = config["samples"].keys()
    )
  output:
    os.path.join(OUTDIR, "concatenated_assembly.fasta")
  threads: 1
  resources:
      mem = "5G",
      walltime = "10:30:00",
      nodetype = config["nodetype"]
  shell:
    """
      cat {input} > {output}
    """
 
rule subset_contig_bin:
  input: os.path.join(OUTDIR,"chunked_assembly/{fragment_size}kb/{sample}-{fragment_size}kb.fasta")
  output: os.path.join(OUTDIR,"chunked_assembly/{fragment_size}kb/{sample}-{fragment_size}.tsv")
  threads: 1
  resources:
      mem = "5G",
      walltime = "00:30:00",
      nodetype = config["nodetype"]
  shell:
      """
      grep ">" {input} | sed 's/>//g' | awk '{{ print $1 "\t" "{wildcards.sample}" }}' > {output}
      """

rule concat_contig_bins:
  input:
    expand(os.path.join(OUTDIR,"chunked_assembly/{fragment_size}kb/{sample}-{fragment_size}.tsv"), fragment_size = config["fragment_size"], sample = config["samples"].keys())
  output:
    os.path.join(OUTDIR, "contig_bin_{fragment_size}kb.tsv"),
  threads: 1
  resources:
      mem = "5G",
      walltime = "00:30:00",
      nodetype = config["nodetype"]
  shell:
    """
    cat {input} > {output}
    """
  
