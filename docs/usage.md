# nf-core/detaxizer: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/detaxizer/usage](https://nf-co.re/detaxizer/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

nf-core/detaxizer is a pipeline to assess raw (meta)genomic data for contaminations and filter reads which were classified as contamination. The default taxon classified as contamination is **_Homo sapiens_**.

## Benchmark

Benchmarking with an [artificial real metagenomic reads dataset](https://doi.org/10.5281/zenodo.10472795) is described in this [publication](https://doi.org/10.1101/2025.03.27.645632) with nf-core/detaxizer 1.1.0. The best performing decontamination was achieved by nf-core/detaxizer with the combination of bbduk with GRCh38 AWS igenome and Kraken2 with the Kraken2 Standard database. This setting reached a recall of 0.99962 (3,770 false negatives of 10,000,002 human read pairs) but a precision of 0.99150 (85,741 false positives of 21,172,961 microbial read pairs). The following command mirrors the settings in the benchmark (but with version 1.2.0 instead of 1.1.0):

```bash
NXF_VER=24.04.4 nextflow run nf-core/detaxizer -r 1.2.0 -profile singularity --input samplesheet.csv --classification_bbduk --classification_kraken2 --outdir results_recall
```

To best retention of microbial reads (precision of 0.99922 = 7,654 false positives of 21 million microbial read pairs) at the cost of higher non-detected human reads (recall of 0.99303 = 69,725 false negatives of 10 million human read pairs) was achieved in the [benchmark](https://doi.org/10.1101/2025.03.27.645632) with the following settings (translated from version 1.1.0 to version 1.2.0):

```bash
NXF_VER=24.04.4 nextflow run nf-core/detaxizer -r 1.2.0 -profile singularity --input samplesheet.csv --classification_kraken2 --kraken2db https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240605.tar.gz --outdir results_precision
```

> [!TIP]
> Remote databases can complicate caching, so that resuming or repeating will re-download large database files. It is recommended instead to pre-download remote databases and refer to local copies.

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 4 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 4 columns to match those defined in the table below. For single-end short reads use the column `short_reads_fastq_1`, for long reads use the column `long_reads_fastq_1`.

A final samplesheet file consisting of both single- and paired-end data may look something like the one below. This is for 5 samples, showing all possible combinations of short and long reads.

```csv title="samplesheet.csv"
sample,short_reads_fastq_1,short_reads_fastq_2,long_reads_fastq_1
SINGLE_END_SHORT,AEG588A1_S1_L002_R1_001.fastq.gz,,
PAIRED_END_SHORT,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz,
SINGLE_END_LONG,,,AEG588A3_001.fastq.gz
SINGLE_END_SHORT_LONG,AEG588A4_S1_L002_R1_001.fastq.gz,,AEG588A4_001.fastq.gz
PAIRED_END_PLUS_LONG,AEG588A5_S1_L002_R1_001.fastq.gz,AEG588A5_S1_L002_R2_001.fastq.gz,AEG588A5_001.fastq.gz
```

| Column                | Description                                                                                                                                                                            |
| --------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`              | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `short_reads_fastq_1` | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or "fq.gz". Optional, if `long_reads_fastq_1` is also provided.          |
| `short_reads_fastq_2` | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or "fq.gz". Optional. Only used for paired-end files.                    |
| `long_reads_fastq_1`  | Full path to FastQ file for long reads. File has to be gzipped and have the extension ".fastq.gz" or "fq.gz". Optional. Use only for long reads.                                       |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Databases

The databases used by detaxizer have an influence on the amount of false positives (classified as contamination although not originating from that taxon/taxa) and false negatives (not classified as the taxon/taxa defined as contamination although of such origin).

The task of decontamination has to be balanced out between false positives and false negatives depending on what is needed in your use case.

> [!NOTE]
> Be aware that the `tax2filter` (default _Homo sapiens_) has to be in the provided kraken2 database (if kraken2 is used) and that the reference for bbduk (provided by the `fasta_bbduk` parameter) should contain the taxa to filter/assess if it is wanted to assess/remove the same taxa as in `tax2filter`. This overlap in the databases is not checked by the pipeline. To filter out/assess taxa with bbduk only, the `tax2filter` parameter is not needed but a fasta file with references of these taxa has to be provided.

> [!TIP]
> The choice of database matters for performance and computational requirements as described in our [benchmarking study](https://doi.org/10.1101/2025.03.27.645632).

> [!TIP]
> Local and remote databases can be used. However, remote databases can complicate caching, so that resuming or repeating will re-download large database files (indicated by `staging foreign file`). It is recommended instead to pre-download remote databases and refer to local copies.

### kraken2

For optimal decontamination performance a large kraken2 database should be used. This comes at costs in terms of hardware requirements and download volumes. The default is a large database called "Standard", with ~60GB size and requires ~80GB RAM. For the largest kraken2 standard database (which can be found [here](https://benlangmead.github.io/aws-indexes/k2)) at least 100 GB of memory should be available, depending on the size of your data the required memory may be higher. For standard decontamination tasks the Standard-8 GB database (i.e. caped at 8GB) can be used, but it should always be kept in mind that this may lead to more false negatives (but conversely potentially less false positives).

To build your own database refer to [this site](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#custom-databases).

### bbduk

bbduk uses a fasta file which contains sequences from the taxon/taxa classified as contamination. Default is the `GRCh38` human reference genome. Provide a custom file using the `fasta_bbduk` parameter.

### blastn

The blastn database is built from a fasta file. Default is the `GRCh38` human reference genome. To decrease the amount of false negatives in this step or include different taxa, a database of several taxa can be used. The fasta containing desired sequences has to be provided to the pipeline by using the `fasta_blastn` parameter.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/detaxizer --input ./samplesheet.csv --outdir ./results -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/detaxizer -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

Before and after the execution of the pipeline the headers inside the `.fastq.gz` files are renamed. This step is necessary to avoid difficulties with different header formats in the pipeline. The renamed headers will never be shown to you, except when looking into the work directory. Only the original read headers are shown in the results.

To change the taxon or taxonomic subtree which is classified by kraken2 as contamination use the `tax2filter` parameter (default `Homo sapiens`). The taxon has to be in the kraken2 database used, which can be specified using the `kraken2db` parameter.

To change what is classified by `bbduk`, a fasta containing the sequences of the contaminant taxon/taxa has to be provided using the `fasta_bbduk` parameter.

If you want to run `bbduk` use the `--classification_bbduk` flag. For running both classification steps and use the merged output for filtering, use both flags (`--classification_kraken2` and `--classification_bbduk`).

To change the organism(s) which should be validated as contamination(s) by blasting against a database, you have to provide a fasta from which the blastn database is built using the `fasta_blastn` parameter. Also, if just one reference genome is needed for blastn and it is in `igenomes.config` use the according name (e.g. `'GRCh38'`) as `genome` parameter.

blastn can be turned on using the `validation_blastn` parameter.

There are two options for the input of the filter, either the raw reads or the preprocessed ones. The first is the default option. Also, for the definition of the reads to be filtered by their IDs two options are available. Either the default is taken, the output from the classification step (kraken2), or using the output from the `blastn` step. The filtering step can be omitted altogether with `--skip_filter`.

If you want to output the removed reads, use `--output_removed_reads`.

Optional classification of the filtered (and removed) reads can be done using `--classification_kraken2_post_filtering`. This uses the kraken2 database provided by `kraken2db`.

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/detaxizer
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/detaxizer releases page](https://github.com/nf-core/detaxizer/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
