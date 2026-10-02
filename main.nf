#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/detaxizer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/detaxizer
    Website: https://nf-co.re/detaxizer
    Slack  : https://nfcore.slack.com/channels/detaxizer
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DETAXIZER              } from './workflows/detaxizer'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_detaxizer_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_detaxizer_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_DETAXIZER {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    DETAXIZER (
        samplesheet,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir,
    )
    emit:
    multiqc_report = DETAXIZER.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params {

    // Path to comma-separated file containing information about the samples in the experiment.
    input: Path?

    // The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.
    outdir: String?

    // Email address for completion summary.
    email: String?

    // MultiQC report title. Printed as page header, used for filename if not otherwise specified.
    multiqc_title: String?

    // If preprocessing with fastp should be turned on.
    preprocessing: Boolean

    // Signifies that bbduk is used in the classification process. Can be combined with the 'classification_kraken2' parameter to run both.
    classification_bbduk: Boolean

    // Signifies that kraken2 is used in the classification process. Can be combined with the 'classification_bbduk' parameter to run both. For kraken2 alone no parameter is needed.
    classification_kraken2: Boolean

    // If a validation of the classified reads via blastn should be carried out.
    validation_blastn: Boolean

    // If the filtered reads should be classified with kraken2.
    classification_kraken2_post_filtering: Boolean

    // When a validation via blastn is wanted but the filtering should use the IDs from the classification process.
    filter_with_classification: Boolean

    // If the filtering step should be skipped.
    skip_filter: Boolean

    // Select the read-filtering tool: seqkit or bbmap. seqkit normalizes FASTQ headers by temporarily renaming them; bbmap uses filterbyname.sh for exact header matching -- Note: BBTools I/O forces any base that is N to Q=0 (!).
    filtering_tool: String = 'seqkit'

    // If the removed reads should also be written to the output folder.
    output_removed_reads: Boolean

    // If the pre-processed reads should be used by the filter.
    filter_trimmed: Boolean

    // Save intermediates to the results folder.
    save_intermediates: Boolean

    // Location of the fasta which contains the contaminant sequences.
    fasta_bbduk: Path?

    // Length of k-mers for classification carried out by bbduk
    bbduk_kmers: Integer = 27

    // The database which is used in the classification step. Please be aware that this default database will require ~60GB download and ~80GB RAM.
    kraken2db: Path = 'https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240605.tar.gz'

    // Save unclassified reads and classified reads (those assigned to any taxon, not specifically assessed or filtered) to separate files.
    save_output_fastqs: Boolean

    // Save unclassified reads and classified reads (those assigned to any taxon, not specifically assessed or filtered) to separate files. For the filtered reads.
    save_output_fastqs_filtered: Boolean

    // Save unclassified reads and classified reads (those assigned to any taxon, not specifically assessed or filtered) to separate files. For the removed reads.
    save_output_fastqs_removed: Boolean

    // Confidence in the classification of a read as a certain taxon.
    kraken2confidence: Float = 0.0

    // Confidence in the classification of a read as a certain taxon. For the filtered reads.
    kraken2confidence_filtered: Float = 0.0

    // Confidence in the classification of a read as a certain taxon. For the removed reads.
    kraken2confidence_removed: Float = 0.0

    // If a read has less k-mers assigned to the taxon/taxa to be assessed/to be filtered the read is ignored by the pipeline.
    cutoff_tax2filter: Integer = 0

    // Ratio per read of assigned to tax2filter k-mers to k-mers assigned to any other taxon (except unclassified).
    cutoff_tax2keep: Float = 0.0

    // Ratio per read of assigned to tax2filter k-mers to unclassified k-mers.
    cutoff_unclassified: Float = 0.0

    // The taxon or taxonomic group to be assessed or filtered by the pipeline.
    tax2filter: String = 'Homo sapiens'

    // Location of the fasta from which the blastn database will be constructed.
    fasta_blastn: Path?

    // Coverage is the percentage of the query sequence which can be found in the alignments of the sequence match. It can be used to fine-tune the validation step.
    blast_coverage: Float = 40.0

    // The expected(e)-value contains information on how many hits of the same score can be found in a database of the size used in the query by chance. The parameter can be used to fine-tune the validation step.
    blast_evalue: Float = 0.01

    // Identity is the percentage of the exact matches in the query and the sequence found in the database. The parameter can be used to fine-tune the validation step.
    blast_identity: Float = 40.0

    // fastp option defining the minimum readlength of a read
    reads_minlength: Integer = 0

    // fastp option defining if the reads which failed to be trimmed should be saved
    fastp_save_trimmed_fail: Boolean

    // fastp option to define the threshold of quality of an individual base
    fastp_qualified_quality: Integer = 0

    // fastp option to define the mean quality for trimming
    fastp_cut_mean_quality: Integer = 1

    // fastp option if duplicates should be filtered or not before classification
    fastp_eval_duplication: Boolean

    // fastp option to define if the clipped reads should be saved
    save_clipped_reads: Boolean

    // Name of iGenomes reference.
    genome: String?

    // Do not load the iGenomes reference config.
    igenomes_ignore: Boolean

    // Save the reference genome and its index files in the results directory.
    saveReference: Boolean

    // The base path to the igenomes reference files
    igenomes_base: String?

    // Turn on generation of samplesheets for downstream pipelines.
    generate_downstream_samplesheets: Boolean

    // Specify a comma separated string in quotes to specify which pipeline to generate a samplesheet for.
    generate_pipeline_samplesheets: String = 'taxprofiler,mag'

    // Git commit id for Institutional configs.
    custom_config_version: String?

    // Base directory for Institutional configs.
    custom_config_base: String?

    // Institutional config name.
    config_profile_name: String?

    // Institutional config description.
    config_profile_description: String?

    // Institutional config contact information.
    config_profile_contact: String?

    // Institutional config URL link.
    config_profile_url: String?

    // Display version and exit.
    version: Boolean

    // Method used to save pipeline results to output directory.
    publish_dir_mode: String?

    // Email address for completion summary, only when pipeline fails.
    email_on_fail: String?

    // Send plain-text email instead of HTML.
    plaintext_email: Boolean

    // File size limit when attaching MultiQC reports to summary emails.
    max_multiqc_email_size: String = '25.MB'

    // Do not use coloured log outputs.
    monochrome_logs: Boolean

    // Custom config file to supply to MultiQC.
    multiqc_config: Path?

    // Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file
    multiqc_logo: Path?

    // Custom MultiQC yaml file containing HTML including a methods description.
    multiqc_methods_description: Path?

    // Boolean whether to validate parameters against the schema at runtime
    validate_params: Boolean = true

    // Base URL or local path to location of pipeline test dataset files
    pipelines_testdata_base_path: String?

    // Suffix to add to the trace report filename. Default is the date and time in the format yyyy-MM-dd_HH-mm-ss.
    trace_report_suffix: String

    // Display the help message.
    help: Boolean

    // Display the full detailed help message.
    help_full: Boolean

    // Display hidden parameters in the help message (only works when --help or --help_full are provided).
    show_hidden: Boolean
}

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_DETAXIZER (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_DETAXIZER.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
