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
    input: Path?
    outdir: String?
    email: String?
    multiqc_title: String?
    preprocessing: Boolean
    classification_bbduk: Boolean
    classification_kraken2: Boolean
    validation_blastn: Boolean
    classification_kraken2_post_filtering: Boolean
    filter_with_classification: Boolean
    skip_filter: Boolean
    filtering_tool: String = 'seqkit'
    output_removed_reads: Boolean
    filter_trimmed: Boolean
    save_intermediates: Boolean
    fasta_bbduk: Path?
    bbduk_kmers: Integer = 27
    kraken2db: Path = 'https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240605.tar.gz'
    save_output_fastqs: Boolean
    save_output_fastqs_filtered: Boolean
    save_output_fastqs_removed: Boolean
    kraken2confidence: Float = 0.0
    kraken2confidence_filtered: Float = 0.0
    kraken2confidence_removed: Float = 0.0
    cutoff_tax2filter: Integer = 0
    cutoff_tax2keep: Float = 0.0
    cutoff_unclassified: Float = 0.0
    tax2filter: String = 'Homo sapiens'
    fasta_blastn: Path?
    blast_coverage: Float = 40.0
    blast_evalue: Float = 0.01
    blast_identity: Float = 40.0
    reads_minlength: Integer = 0
    fastp_save_trimmed_fail: Boolean
    fastp_qualified_quality: Integer = 0
    fastp_cut_mean_quality: Integer = 1
    fastp_eval_duplication: Boolean
    save_clipped_reads: Boolean
    genome: String?
    igenomes_ignore: Boolean
    saveReference: Boolean
    igenomes_base: String?
    generate_downstream_samplesheets: Boolean
    generate_pipeline_samplesheets: String = 'taxprofiler,mag'
    custom_config_version: String?
    custom_config_base: String?
    config_profile_name: String?
    config_profile_description: String?
    config_profile_contact: String?
    config_profile_url: String?
    version: Boolean
    publish_dir_mode: String?
    email_on_fail: String?
    plaintext_email: Boolean
    max_multiqc_email_size: String = '25.MB'
    monochrome_logs: Boolean
    multiqc_config: Path?
    multiqc_logo: Path?
    multiqc_methods_description: Path?
    validate_params: Boolean = true
    pipelines_testdata_base_path: String?
    trace_report_suffix: String
    help: Boolean
    help_full: Boolean
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
