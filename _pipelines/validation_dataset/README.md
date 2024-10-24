# Validation Datasets
This directory contains the scripts to create the validation datasets from the isolate data + ZymoHMW. The validation datasets are used to evaluate the performance of SemiBin3.

## Scripts
Before running the snakemake pipelines the `find_samples.ipynb` script was used to create a `samples.yaml` file that contains the sample names and the path to the fastq and assembly file.

Furthermore, the `find_samples.ipynb` can create a subsampling schema, where the isolates are subsampled to a certain coverage following a log-normal distribution. The subsampling schema is used to create another type of validation dataset.

### Create validation datasets
Current ways of constructing validation datasets are:
1. Simple concatenation of isolates and ZymoHMW
2. Subsampling isolates to a certain coverage following a log-normal distribution

#### workflows/concatenate_isolates.smk
This workflow simply concatenates all the assemblies and fastq files and then creates a new `pileup.bed` file. With these files the `SemiBin3` workflow can be run.

#### workflows/sample_qc & workflows/subsample.smk
The `sample_qc` workflow filters the fastq files to remove reads with a mean quality score below `q=10` and reads with a length below `200bp`. The `subsample.smk` workflow subsamples the isolates to a certain coverage following a log-normal distribution. After subsampling the fastq-file is used for assembly with `flye` and a new `pileup.bed` file is created.


#### Creating the pileup.bed file
The `pileup.bed` file is created by running the `pileup.smk` workflow.



