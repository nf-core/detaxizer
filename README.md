<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-detaxizer_logo_dark.png">
    <img alt="nf-core/detaxizer" src="docs/images/nf-core-detaxizer_logo_light.png">
  </picture>
</h1>

[![Cite Publication](https://img.shields.io/badge/Cite%20Us!-Cite%20Publication-important?labelColor=000000)](https://doi.org/10.1093/nargab/lqaf125)

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://github.com/codespaces/new/nf-core/detaxizer)
[![GitHub Actions CI Status](https://github.com/nf-core/detaxizer/actions/workflows/nf-test.yml/badge.svg)](https://github.com/nf-core/detaxizer/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/detaxizer/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/detaxizer/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/detaxizer/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.10877147-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.10877147)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.8-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.4.1-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.4.1)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/detaxizer)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23detaxizer-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/detaxizer)[![Follow on Bluesky](https://img.shields.io/badge/bluesky-%40nf__core-1185fe?labelColor=000000&logo=bluesky)](https://bsky.app/profile/nf-co.re)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/detaxizer** is a bioinformatics pipeline that checks for the presence of a specific taxon in (meta)genomic fastq files and to filter out this taxon or taxonomic subtree. The process begins with quality assessment via FastQC and optional preprocessing (adapter trimming, quality cutting and optional length and quality filtering) using fastp, followed by taxonomic classification with kraken2 and/or bbduk, and optionally employs blastn for validation of the reads associated with the identified taxa. Users must provide a samplesheet to indicate the fastq files and, if utilizing bbduk in the classification and/or the validation step, fasta files for usage of bbduk and creating the blastn database to verify the targeted taxon.

![detaxizer metro workflow](docs/images/Detaxizer_metro_workflow.png)

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Optional pre-processing ([`fastp`](https://github.com/OpenGene/fastp))
3. Classification of reads ([`Kraken2`](https://ccb.jhu.edu/software/kraken2/), and/or [`bbduk`](https://sourceforge.net/projects/bbmap/))
4. Optional validation of searched taxon/taxa ([`blastn`](https://blast.ncbi.nlm.nih.gov/Blast.cgi))
5. Filtering of the searched taxon/taxa from the reads (either from the raw files or the preprocessed reads, using either the output from the classification (kraken2 and/or bbduk) or blastn)
6. Summary of the processes (how many were classified and optionally how many were validated)
7. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

```csv title="samplesheet.csv"
sample,short_reads_fastq_1,short_reads_fastq_2,long_reads_fastq_1
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,AEG588A1_S1_L002_R3_001.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end). A third fastq file can be provided if long reads are present in your project. For more detailed information about the samplesheet, see the [usage documentation](docs/usage.md).

Now, you can run the pipeline using:

```bash
nextflow run nf-core/detaxizer \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --classification_bbduk \
   --classification_kraken2 \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/detaxizer/usage) and the [parameter documentation](https://nf-co.re/detaxizer/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/detaxizer/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/detaxizer/output).

Generated samplesheets from the directory `/downstream_samplesheets/` can be used for the pipelines:

- [nf-core/mag](https://nf-co.re/mag)
- [nf-core/taxprofiler](https://nf-co.re/taxprofiler)

## Credits

nf-core/detaxizer was originally written by [Jannik Seidel](https://github.com/jannikseidelQBiC) at the [Quantitative Biology Center (QBiC)](http://qbic.life/).

We thank the following people for their extensive assistance in the development of this pipeline:

- [Daniel Straub](https://github.com/d4straub)

This work was initially funded by the German Center for Infection Research (DZIF).

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#detaxizer` channel](https://nfcore.slack.com/channels/detaxizer) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use nf-core/detaxizer for your analysis, please cite it using the following publication:

> **nf-core/detaxizer: a benchmarking study for decontamination from human sequences.**
>
> Jannik Seidel, Camill Kaipf, Daniel Straub, Sven Nahnsen
>
> _NAR Genom and Bioinform._ 2025;7(3):lqaf125, doi: [10.1093/nargab/lqaf125](https://doi.org/10.1093/nargab/lqaf125).

Additionally, the following doi can be cited: [10.5281/zenodo.10877147](https://doi.org/10.5281/zenodo.10877147)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
