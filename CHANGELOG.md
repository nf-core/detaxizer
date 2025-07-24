# nf-core/detaxizer: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.2.0dev - [DATE]

### `Added`

- [PR #70](https://github.com/nf-core/detaxizer/pull/70) - Filtering is now default, `--skip_filter` was added
- [PR #71](https://github.com/nf-core/detaxizer/pull/71) - Add usage information learned from our benchmarking

### `Changed`

- [PR #65](https://github.com/nf-core/detaxizer/pull/65),[PR #69](https://github.com/nf-core/detaxizer/pull/69) - Template update for nf-core/tools 3.3.2 (by @d4straub)

### `Fixed`

### `Dependencies`

| Software | Previous version | New version |
| -------- | ---------------- | ----------- |
| MultiQC  | 1.27             | 1.29        |

### `Deprecated`

- [PR #70](https://github.com/nf-core/detaxizer/pull/70) - Filtering is now default, `--enable_filter` was removed and replaced by `--skip_filter`

## v1.1.0 - Kombjuudr - [2024-11-08]

### `Added`

- [PR #34](https://github.com/nf-core/detaxizer/pull/34) - Added bbduk to the classification step (kraken2 as default, both can be run together) (by @jannikseidelQBiC)
- [PR #34](https://github.com/nf-core/detaxizer/pull/34) - Added `--fasta_bbduk` parameter to provide a fasta file with contaminants (by @jannikseidelQBiC)
- [PR #34](https://github.com/nf-core/detaxizer/pull/34) - Rewrote summary step of classification to be usable with bbduk and/or kraken2 (by @jannikseidelQBiC)
- [PR #34](https://github.com/nf-core/detaxizer/pull/34) - Made preprocessing with fastp optional and added the parameter `--fastp_eval_duplication` to turn on duplication removal (off as default, was on/not changeable in v1.0.0) (by @jannikseidelQBiC)
- [PR #34](https://github.com/nf-core/detaxizer/pull/34) - Optionally the removed reads can now be written to the output folder (by @jannikseidelQBiC)
- [PR #34](https://github.com/nf-core/detaxizer/pull/34) - Added optional classification of filtered and removed reads via kraken2 (by @jannikseidelQBiC)
- [PR #39](https://github.com/nf-core/detaxizer/pull/39) - Added generation of input samplesheet for nf-core/mag, nf-core/taxprofiler (by @Joon-Klaps)

#### Parameters

Added parameters:

| Parameter                                 |
| ----------------------------------------- |
| `--fasta_bbduk`                           |
| `--preprocessing`                         |
| `--output_removed_reads`                  |
| `--classification_kraken2`                |
| `--classification_bbduk`                  |
| `--kraken2confidence_filtered`            |
| `--kraken2confidence_removed`             |
| `--classification_kraken2_post_filtering` |
| `--fastp_eval_duplication`                |
| `--bbduk_kmers`                           |

Changed default values of parameters:

| Parameter                  | Old default value                                                             | New default value                                                             |
| -------------------------- | ----------------------------------------------------------------------------- | ----------------------------------------------------------------------------- |
| `--fastp_cut_mean_quality` | 15                                                                            | 1                                                                             |
| `--kraken2db`              | 'https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20231009.tar.gz' | 'https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240605.tar.gz' |
| `--kraken2confidence`      | 0.05                                                                          | 0.00                                                                          |
| `--tax2filter`             | 'Homo'                                                                        | 'Homo sapiens'                                                                |
| `--cutoff_tax2filter`      | 2                                                                             | 0                                                                             |
| `--cutoff_tax2keep`        | 0.5                                                                           | 0.0                                                                           |

### `Changed`

- [PR #42](https://github.com/nf-core/detaxizer/pull/42) - Template update for nf-core/tools 3.0.2, for details read [this blog post](https://nf-co.re/blog/2024/tools-3_0_0#important-template-updates)

### `Fixed`

- [PR #33](https://github.com/nf-core/detaxizer/pull/33) - Addition of quotation marks in `parse_kraken2report.nf` prevents failure of the pipeline when using a taxon with space (e.g. Homo sapiens) with the `--tax2filter` parameter (by @jannikseidelQBiC)
- [PR #34](https://github.com/nf-core/detaxizer/pull/34) - Made validation via blastn optional by default (by @jannikseidelQBiC)
- [PR #34](https://github.com/nf-core/detaxizer/pull/34) - Changed parameter `--fasta` to `--fasta_blastn` (by @jannikseidelQBiC)

### `Dependencies`

Updated and added dependencies

| Tool    | Previous version | Current version |
| ------- | ---------------- | --------------- |
| bbmap   | -                | 39.10           |
| blastn  | 2.14.1           | 2.15.0          |
| multiQC | 1.21             | 1.25.1          |
| kraken2 | 2.1.2            | 2.1.3           |
| seqkit  | 2.8.0            | 2.8.2           |

### `Deprecated`

| Parameter       | New parameter         | Reason                                                                                                                                              |
| --------------- | --------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--fasta`       | `--fasta_blastn`      | Introduction of fasta_bbduk; necessary to further distinguish the two parameters                                                                    |
| `--skip_blastn` | `--validation_blastn` | blastn is now to be enabled on purpose; too resource intensive for a default setting                                                                |
| `--max_cpus`    | -                     | New behavior of [nextflow](https://www.nextflow.io/docs/latest/reference/process.html#resourcelimits), `resourceLimits` can now be set via a config |
| `--max_memory`  | -                     | New behavior of [nextflow](https://www.nextflow.io/docs/latest/reference/process.html#resourcelimits), `resourceLimits` can now be set via a config |
| `--max_time`    | -                     | New behavior of [nextflow](https://www.nextflow.io/docs/latest/reference/process.html#resourcelimits), `resourceLimits` can now be set via a config |

## v1.0.0 - Kobbfarbad - [2024-03-26]

Initial release of nf-core/detaxizer, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Added preprocessing part of workflow (tool: fastp)
- Added examination part of workflow (tools: kraken2 and optionally blastn) and a summary step for this part
- Added optional filtering step (either applied to raw reads or preprocessed ones, and either with the output of kraken2 or blastn)
- Added read renaming step at the beginning and end of the workflow to cope with different fastq header formats while keeping original headers in the filtered fastqs

### `Fixed`

### `Dependencies`

### `Deprecated`
