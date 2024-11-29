
import os
import glob
import datetime
import pandas as pd

SNAKEDIR = os.path.dirname(workflow.snakefile)
configfile: os.path.join(SNAKEDIR, "..", "config", "config-fragmentation_benchmark.yaml")

OUTDIR = config["outdir"]
SAMPLES = config["samples"].keys()


rule all:
  input:
    os.path.join(OUTDIR, "concatenated_assembly.fasta")

################################################################################
# Monoculture data subsetting
# rule extract_reads:
#   input:
#     lambda wildcards: expand(config["samples"][wildcards.sample]["coverage"], sample = config["samples"].keys()),
#   output:
#     temp(os.path.join(OUTDIR, "fastq_files", "{sample}.fastq"))
#   threads: 1
#   resources:
#       mem = "5G",
#       walltime = "10:30:00",
#       nodetype = config["nodetype"]
#   shell:
#     """
#       samtools fastq -T MM,ML {input} > {output}
#     """
  
  
 
# rule concatenate_reads:
#   input:
#     i = expand(os.path.join(OUTDIR, "fastq_files", "{sample}.fastq"), sample = SAMPLES),
#   output:
#     os.path.join(OUTDIR, "concatenated_reads.fastq"),
#   threads: 1
#   resources:
#       mem = "5G",
#       walltime = "10:30:00",
#       nodetype = config["nodetype"]
#   shell:
#     """
#       cat {input} > {output}
#     """

rule rename_assemblies:
  input:
    lambda wildcards: config["samples"][wildcards.sample]["assembly"]
  output:
    temp(os.path.join(OUTDIR, "assembly_files", "{sample}.fasta"))
  threads: 1
  resources:
      mem = "5G",
      walltime = "10:30:00",
      nodetype = config["nodetype"]
  shell:
    """
      cp {input} {output}
      sed '/^>/ s/^>/>{wildcards.sample}_/' {output} > {output}.tmp
      mv {output}.tmp {output}
      
    """


rule concatenate_assembly:
  input:
    expand(
      os.path.join(OUTDIR, "assembly_files", "{sample}.fasta"),
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
 


      
# rule map_mod_calls:
#     conda: 
#         "envs/mapreads.yaml"
#     input:
#         assembly = os.path.join(OUTDIR, "concatenated_assembly.fasta"),
#         reads = os.path.join(OUTDIR, "concatenated_reads.fastq"),
#     output:
#         bam = os.path.join(OUTDIR, "mod.bam")
#     threads: 40,
#     params:
#         mem = "50g",
#         ref = "10g"
#     resources:
#         mem = "100G",
#         nodetype="high-mem",
#         walltime="2-00:00:00",
#     shell:
#         """
#         minimap2 -ax map-ont -t {threads} -I {params.mem} -K {params.ref} -y {input.assembly} {input.reads} | \
#         samtools view -bS | \
#         samtools sort > {output.bam}
#         """


# rule index_bam:
#     conda: 
#         "envs/mapreads.yaml"
#     input:
#         bam = os.path.join(OUTDIR, "mod.bam")
#     output:
#         bai = os.path.join(OUTDIR, "mod.bam.bai")
#     threads: 10,
#     resources:
#         mem ="10G",
#         nodetype=config["nodetype"],
#         walltime="1-00:00:00",
#     shell:
#         """
#         samtools index -@ {threads} {input.bam}
#         """

# rule genomecov:
#     conda: "envs/bedtools.yaml"
#     input:
#       bam = os.path.join(OUTDIR, "mod.bam")
#     output:
#       txt = os.path.join(OUTDIR, "genomecov.txt")
#     threads: 1,
#     resources:
#         mem ="5G",
#         nodetype=config["nodetype"],
#         walltime="1-00:00:00",
#     shell:
#       """
#         bedtools genomecov -bga -ibam {input} > {output}
#       """
      
