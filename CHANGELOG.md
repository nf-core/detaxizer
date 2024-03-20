# nf-core/detaxizer: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0 - Kobbfarbad - [2024-03-20]

Initial release of nf-core/detaxizer, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Added preprocessing part of workflow (tool: fastp)
- Added examination part of workflow (tools: kraken2 and optionally blastn) and a summary step for this part
- Added optional filtering step (either applied to raw reads or preprocessed ones, and either with the output of kraken2 or blastn)
- Added read renaming step at the beginning and end of the workflow to cope with different fastq header formats while keeping original headers in the filtered fastqs

### `Fixed`

### `Dependencies`

### `Deprecated`
