# nf-core/detaxizer: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0dev - [date]

### `Added`

- Added bbduk to the classification step (kraken2 as default, both can be run together)
- Added `fasta_bbduk` parameter to provide a fasta file with contaminants
- Rewrote summary step of classification to be usable with bbduk and/or kraken2

### `Fixed`

- [PR #33](https://github.com/nf-core/detaxizer/pull/33) - Addition of quotation marks in `parse_kraken2report.nf` prevents failure of the pipeline when using a taxon with space (e.g. Homo sapiens) with the `tax2filter` parameter.
- Made validation via blastn optional by default
- Changed parameter `fasta` to `fasta_blastn`

### `Dependencies`

### `Deprecated`

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
