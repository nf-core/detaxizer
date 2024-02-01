# nf-core/detaxizer: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory. `<sample>` is a placeholder for the real sample name provided in the `samplesheet.csv`.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Renaming](#renaming) - Renaming the read headers before and after the pipeline run to avoid difficulties with different header formats
- [FastQC](#fastqc) - Raw read QC - Output not in the results directory by default
- [fastp](#fastp) - Preprocessing of raw reads
- [kraken2](#kraken2) - Classification of the preprocessed reads and extracting the searched taxa from the results
- [blastn](#blastn) - Validation of the reads classified as the searched taxa and extracting ids of validated reads
- [filter](#filter) - (Optional) filtering of the raw or preprocessed reads using either the read ids from kraken2 output or blastn output
- [summary](#summary) - The summary of the classification and the optional validation
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Renaming

Before and after (if using the filter) the execution of the pipeline the headers inside the `.fastq.gz` files are renamed. As stated above, this step is necessary to avoid difficulties with different header formats in the pipeline. The renamed headers will never be shown to you, except when looking into the work directory. Only the original read headers are shown in the results.

### FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/). The output is not shown in the results folder by default.

### fastp

fastp performs preprocessing of the reads (adapter/quality trimming). For details of the output, please refer to [this site](https://nf-co.re/modules/fastp).

<details markdown="1">
<summary>Output files</summary>

- `fastp/`: Contains the output from the preprocessing step.
  - `<sample>_longReads/`: If long reads are present in your `samplesheet.csv` this folder is generated containing the fastp-report.
    - `<sample>_longReads.fastp.html`: The report on the preprocessing step.
    - `<sample>_longReads.fastp.json`: The data on the preprocessing step in `.json`-format.
  - `<sample>_R1/`: If single-end short reads are present in your `samplesheet.csv` this folder is generated.
    - same pattern as in `<sample>_longReads/` with the prefix `<sample>_R1.fastp.*`.
  - `<sample>/`: For paired-end short reads in your `samplesheet.csv` this folder is generated.
    - same pattern as in `<sample>_longReads/` with the prefix `<sample>.fastp.*`.

</details>

### kraken2

kraken2 classifies the reads. The important files are `*.classifiedreads.txt`, `*.kraken2.report.txt`, `isolated/*.classified.txt` and `summary/*.kraken2_summary.tsv`. The first contains all reads, their classification and how many k-mers were assigned to which taxon. The second contains statistics on how many reads were classified as which taxon. Next is a file which is similar to the first one but only contains the read ids which were classified as the taxon/taxa to assess/to filter together with the whole information from the first file for the individual read ids. Last the summary gives you a fast overview of how many reads were passed to kraken2 and how many were classified as the taxon/taxa to assess or to filter. `<sample>` can be replaced by `<sample>_longReads`, `<sample>_R1` or left as `<sample>` depending on the cases mentioned in [fastp](#fastp).

<details markdown="1">
<summary>Output files</summary>

- `kraken2/`: Contains the output from the classification step.
  - `isolated/`: Contains the isolated lines and ids for the taxon/taxa mentioned in the `tax2filter` parameter.
    - `<sample>.classified.txt`: The whole kraken2 output for the taxon/taxa mentioned in the `tax2filter` parameter.
    - `<sample>.ids.txt`: The ids from the whole kraken2 output assigned to the taxon/taxa mentioned in the `tax2filter` parameter.
  - `summary/`: Summary of the kraken2 process.
    - `<sample>.kraken2_summary.tsv`: Contains two three columns, column 1 is the sample name, column 2 the amount of lines in the untouched kraken2 output and column 3 the amount of lines in the isolated output.
  - `taxonomy/`: Contains the list of taxa to filter/to assess for.
    - `taxa_to_filter.txt`: Contains the taxon ids of all taxa to assess the data for or to filter out.
  - `<sample>.classifiedreads.txt`: The whole kraken2 output for all reads.
  - `<sample>.kraken2.report.txt`: Statistics on how many reads where assigned to which taxon/taxonomic group.

</details>

### blastn

blastn can validate the reads classified by kraken2 as the taxon/taxa to be assessed/to be filtered. To reduce computational burden only the highest scoring hit per input sequence is returned. If in any case one would need more information this can be done via the `max_hsps`- and `max_target_seqs`-flags in the `modules.config`.

<details markdown="1">
<summary>Output files</summary>

- `blast/`
  - `filteredIdentCov/`: The read ids and statistics of the reads which were validated by blastn to be the taxon/taxa to assess/to filter.
    - `<sample>_R1.identcov.txt`: File is present for single-end and paired-end short reads.
    - `<sample>_R2.identcov.txt`: File is present for paired-end short reads.
    - `<sample>_longReads.identcov.txt`: File is present for long reads.
  - `summary/`: Short overview of the amount of reads which were validated by blastn.
    - `<sample>.blastn_summary.tsv`: `<sample>` can be one of two options for this file. Either stay as `<sample>` or be `<sample>_longReads` for long reads.

</details>

### filter

In this folder, the filtered and re-renamed reads can be found. This result has to be carefully examined using the other information in the results folder.

<details markdown="1">
<summary>Output files</summary>

- `filter/`: Folder containing the filtered and re-renamed reads.
  - `<sample>_filtered.fastq.gz`: The filtered reads, `<sample>` can stay as `<sample>` for single-end short reads, take the pattern `<sample>_{R1,R2}` for paired-end reads and `<sample>_longReads` for long reads.

</details>

### summary

The summary file lists all statistics of kraken2 and blastn per sample. It is a combination of the summary files of kraken2 and blastn and can be used for a quick overview of the pipeline run. If blastn is skipped, then only the statistics of kraken2 is shown.

|                                                                                                                    | kraken2                    | isolatedkraken2                         | blastn_unique_ids                                                         | blastn_lines                         | filteredblastn_unique_ids                                                                                                    | filteredblastn_lines                                                               |
| ------------------------------------------------------------------------------------------------------------------ | -------------------------- | --------------------------------------- | ------------------------------------------------------------------------- | ------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------- |
| `<sample>` (For short reads it is the same as in the `samplesheet.csv`, for long reads it is `<sample>_longReads`) | Read IDs in kraken2 output | Read IDs in the isolated kraken2 output | Number of unique IDs in blastn output, should be the same as blastn_lines | Number of lines in the blastn output | Number of IDs in the blastn output after the filtering for identity and coverage, should be the same as filteredblastn_lines | Number of lines in the blastn output after the filtering for identity and coverage |

<details markdown="1">
<summary>Output files</summary>

- `summary/`: Folder containing the summary.
  - `summary.tsv`: File containing the summary in the format stated above.

</details>

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
