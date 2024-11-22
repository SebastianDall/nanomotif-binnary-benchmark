import os
import glob
import datetime

SNAKEDIR = os.path.dirname(workflow.snakefile)
configfile: os.path.join(SNAKEDIR, "..", "config", "config.yaml")

# pod5 view /projects/dark_science/datasets/v1.0.1/pod5/PaPr00000216MP.pod5 --ids --no-header --output all_ids.txt
# head --lines 10000 all_ids.txt > 10k_ids.txt
# pod5 filter /projects/dark_science/datasets/v1.0.1/pod5/PaPr00000216MP.pod5 --ids 10k_ids.txt --output 10k.pod5

rule all:
  input:
    expand(os.path.join("data", "real_communities", "{sample}", "mod.calls.bam"), sample = config["samples"].keys()),

rule basecall_pod5:
    input: lambda wildcards: config["samples"][wildcards.sample]
    output: os.path.join("data", "real_communities", "{sample}", "mod.calls.bam"),
    threads: 5
    resources:
        mem="25G",
        nodetype="gpu",
        gpu_config=config["gpu_config"],
        walltime="2-00:00:00"
    params:
      DORADO=os.path.join(SNAKEDIR, "src", "dorado-0.8.3-linux-x64", "bin", "dorado"),
      # MOD_MODELS=f"{os.path.join(SNAKEDIR, "src", "models", "dna_r10.4.1_e8.2_400bps_sup@v5.0.0_4mC_5mC@v2")},{os.path.join(SNAKEDIR, "src", "models", "dna_r10.4.1_e8.2_400bps_sup@v5.0.0_6mA@v2")}",
      MOD_MODELS = ",".join([
          os.path.join(SNAKEDIR, "src", "models", "dna_r10.4.1_e8.2_400bps_sup@v5.0.0_4mC_5mC@v2"),
          os.path.join(SNAKEDIR, "src", "models", "dna_r10.4.1_e8.2_400bps_sup@v5.0.0_6mA@v2")
      ]),
      MODEL=os.path.join(SNAKEDIR,"src","models", "dna_r10.4.1_e8.2_400bps_sup@v5.0.0"),
      KIT="SQK-NBD114-24"
    shell:
      """
        {params.DORADO} basecaller {params.MODEL} {input} \
          --no-trim \
          --verbose \
          --kit-name {params.KIT} \
          --device cuda:0 > {output}
      """
          # --modified-bases-models {params.MOD_MODELS} \
