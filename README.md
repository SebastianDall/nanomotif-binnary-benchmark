# Nanomotif development
Additional code for reproducing results in "Nanomotif: Leveraging DNA Methylation Motifs for Genome Recovery and Host Association of Plasmids in Metagenomes from Complex Microbial Communities".

## Installation
```bash
mamba create -n nanomotif_dev python=3.12
mamba install -c bioconda nanomotif=0.5.0
```

> This repo assumes understanding of running snakemake pipelines. All pipelines were run with snakemake 7.32. All pipelines are configured using config files. There should be changed to match the paths you are using.

## Creating benchmark dataset
The benchmark dataset was created by running the snakemake pipeline found at `_pipelines/validation_dataset/workflows/fragmentation_benchmark.smk`. This will create a `concatenated_assembly.fasta` and `concatenated_reads.fastq` from the monocultures specified at: `_pipelines/validation_dataset/config/config-fragmentation_benchmark.yaml`.

Each contig in the monoculture assemblies have been fragmented into 20kbp fragments during the creating of `concatenated_assembly.fasta`. To combine fragments to a 1600 kbp fragment run the following:

```bash
python src/concatenate_windows.py <assembly_file> <out_assembly>
```

This will create an assembly with one 1600kbp fragment for each monoculture and several 20kbp fragments.

### Creating pileup
After creating the fragmented assembly, a pileup can be generated using the `_pipelines/validation_dataset/workflows/pileup.smk`.


## Running simulated benchmark
The `_pipelines/simulated_benchmark/workflows/contaminations.smk` and `_pipelines/simulated_benchmark/workflows/include.smk` will create benchmark files where contigs are either shuffled (contamination) or removed (inclusion).

The output will be at `output/benchmark/<detect_contamination|include_contigs>/<benchmark>`. To reproduce the figures run the Rscripts `src/simulated_contamination_benchmark.R` and `src/simulated_include_benchmark.R`.

> Remember to change the paths at the top of the script.




