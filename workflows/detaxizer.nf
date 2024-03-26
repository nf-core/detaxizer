/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_detaxizer_pipeline'
include { getGenomeAttribute     } from '../subworkflows/local/utils_nfcore_detaxizer_pipeline'

include { FASTP             } from '../modules/nf-core/fastp/main'
include { KRAKEN2_KRAKEN2   } from '../modules/nf-core/kraken2/kraken2/main'
include { BLAST_BLASTN      } from '../modules/nf-core/blast/blastn/main'
include { BLAST_MAKEBLASTDB } from '../modules/nf-core/blast/makeblastdb/main'

include { RENAME_FASTQ_HEADERS_PRE              } from '../modules/local/rename_fastq_headers_pre'
include { KRAKEN2PREPARATION                    } from '../modules/local/kraken2preparation'
include { PARSE_KRAKEN2REPORT                   } from '../modules/local/parse_kraken2report'
include { ISOLATE_IDS_FROM_KRAKEN2_TO_BLASTN    } from '../modules/local/isolate_ids_from_kraken2_to_blastn'
include { PREPARE_FASTA4BLASTN                  } from '../modules/local/prepare_fasta4blastn'
include { FILTER_BLASTN_IDENTCOV                } from '../modules/local/filter_blastn_identcov'
include { FILTER                                } from '../modules/local/filter'
include { RENAME_FASTQ_HEADERS_AFTER            } from '../modules/local/rename_fastq_headers_after'
include { SUMMARY_KRAKEN2                       } from '../modules/local/summary_kraken2'
include { SUMMARY_BLASTN                        } from '../modules/local/summary_blastn'
include { SUMMARIZER                            } from '../modules/local/summarizer'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// speficy the fasta channel if it is not provided via --fasta
def fasta = Channel.empty()

if (!params.fasta && !params.skip_blastn) {
    fasta = Channel.fromPath(getGenomeAttribute('fasta'))
} else if (!params.skip_blastn){
    // If params.fasta is there, use it for the creation of the blastn database
    fasta = Channel.fromPath(params.fasta)
}

workflow DETAXIZER {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_samplesheet.branch {
        shortReads: it[1]
        }.set {
            ch_short
        }

    ch_short.shortReads.map{
        meta, short_reads_fastq_1, short_reads_fastq_2, long_reads_fastq_1 ->
            if (short_reads_fastq_2){
                return [meta + [ single_end: false, long_reads: false , amount_of_files: 2 ], [ short_reads_fastq_1, short_reads_fastq_2 ] ]
            } else {
                return [meta + [ id: "${meta.id}_R1", single_end: true, long_reads: false, amount_of_files: 1 ], short_reads_fastq_1 ]
            }
    }.set{
        ch_short
    }

    ch_samplesheet.branch {
        longReads: it[3]
    }.set {
        ch_long
    }

    ch_long.longReads.map {
        meta, short_reads_fastq_1, short_reads_fastq_2, long_reads_fastq_1 ->
            return [meta + [ id: "${meta.id}_longReads", single_end: true, long_reads: true, amount_of_files: 1 ], long_reads_fastq_1 ]
    }.set {
        ch_long
    }

    ch_short_long = ch_short.mix(ch_long)


    //
    // MODULE: Rename Fastq headers
    //
    RENAME_FASTQ_HEADERS_PRE(ch_short_long)


    //
    // MODULE: Run FastQC
    //
    FASTQC (
        RENAME_FASTQ_HEADERS_PRE.out.fastq
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Run fastp
    //
    FASTP (
        RENAME_FASTQ_HEADERS_PRE.out.fastq,
        [],
        params.fastp_save_trimmed_fail,
        []
    )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    //
    // MODULE: Prepare Kraken2 Database
    //
    ch_kraken2_db = Channel.fromPath(params.kraken2db).map {
            item -> [['id': "kraken2_db"], item]
        }
    KRAKEN2PREPARATION (
        ch_kraken2_db
    )
    ch_versions = ch_versions.mix(KRAKEN2PREPARATION.out.versions)

    //
    // MODULE: Run Kraken2
    //

    KRAKEN2_KRAKEN2 (
        FASTP.out.reads,
        KRAKEN2PREPARATION.out.db.first(),
        params.save_output_fastqs,
        true
    )
    ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions.first())

    //
    // MODULE: Parse the taxonomy from the kraken2 report and return all subclasses of the tax2filter
    //
    PARSE_KRAKEN2REPORT(
        KRAKEN2_KRAKEN2.out.report.take(1)
    )
    ch_versions = ch_versions.mix(PARSE_KRAKEN2REPORT.out.versions)

    //
    // MODULE: Isolate the hits for a certain taxa and subclasses
    //
    ch_parsed_kraken2_report = PARSE_KRAKEN2REPORT.out.to_filter.map {meta, path -> path}

    KRAKEN2_KRAKEN2.out.classified_reads_assignment.combine(ch_parsed_kraken2_report).set{ ch_combined }

    ISOLATE_IDS_FROM_KRAKEN2_TO_BLASTN (
        ch_combined
    )

    ch_versions = ch_versions.mix(ISOLATE_IDS_FROM_KRAKEN2_TO_BLASTN.out.versions.first())

    //
    // MODULE: Summarize the kraken2 results and the isolated kraken2 hits
    //
    ch_prepare_summary_kraken2 = KRAKEN2_KRAKEN2.out.classified_reads_assignment.join(ISOLATE_IDS_FROM_KRAKEN2_TO_BLASTN.out.classified).map {
        meta, path1, path2 ->
            return [ meta, [ path1, path2 ] ]
    }

    ch_combined_kraken2 = ch_prepare_summary_kraken2.map {
        meta, path ->
            return [ meta +[ id: meta.id.replaceAll("(_R1|_R2)", "") ] , path]
        }
        .map {
            meta, path ->
                path = path.flatten()
                return [meta, path]
            }

    ch_kraken2_summary = SUMMARY_KRAKEN2(
        ch_combined_kraken2
        )

    ch_versions = ch_versions.mix(ch_kraken2_summary.versions.first())

    // Drop meta of kraken2_summary as it is not needed for the combination step of summarizer
    ch_kraken2_summary = ch_kraken2_summary.summary.map {
            meta, path -> [path]
        }

    if (!params.skip_blastn) {

        //
        // MODULE: Extract the hits to fasta format
        //
        ch_combined = FASTP.out.reads
        .join(
            ISOLATE_IDS_FROM_KRAKEN2_TO_BLASTN.out.classified_ids, by: [0]
        )


        PREPARE_FASTA4BLASTN (
            ch_combined
        )

        ch_versions = ch_versions.mix(PREPARE_FASTA4BLASTN.out.versions.first())

        //
        // MODULE: Run BLASTN
        //
        ch_reference_fasta = fasta

        ch_reference_fasta_with_meta = ch_reference_fasta.map {
            item -> [['id': "id-fasta-for-makeblastdb"], item]
            }

        BLAST_MAKEBLASTDB (
                ch_reference_fasta_with_meta
        )
        ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)

        ch_fasta4blastn = PREPARE_FASTA4BLASTN.out.fasta
            .flatMap { meta, fastaList ->
                if (fastaList.size() == 2) {
                return [
                    [ meta + [ id: "${meta.id}_R1" ], fastaList[0] ],
                    [ meta + [ id: "${meta.id}_R2" ], fastaList[1] ]
                ]

                } else {
                    return [
                        [ meta , fastaList ] ]
                }

            }

        ch_blastn_db = BLAST_MAKEBLASTDB.out.db.first()

        BLAST_BLASTN (
            ch_fasta4blastn,
            ch_blastn_db
        )

        ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first())

        ch_combined_blast = BLAST_BLASTN.out.txt.map {
            meta, path ->
                return [ meta + [ id: meta.id.replaceAll("(_R1|_R2)", "") ], path ]
        }
        .map{
            meta, path -> tuple(groupKey(meta, meta.amount_of_files), path)
        }
        .groupTuple(
                by: [0]
            ).map {
                meta, paths -> [ meta, paths.flatten() ]
                }

        FILTER_BLASTN_IDENTCOV (
            BLAST_BLASTN.out.txt
        )
        ch_versions = ch_versions.mix(FILTER_BLASTN_IDENTCOV.out.versions.first())

        ch_filtered_combined = FILTER_BLASTN_IDENTCOV.out.classified.map {
            meta, path ->
                return [ meta + [ id: meta.id.replaceAll("(_R1|_R2)", "") ], path ]
        }
        .map{
            meta, path -> tuple(groupKey(meta, meta.amount_of_files), path)
        }
        .groupTuple (by: [0])
        .map {
            meta, paths ->
                paths = paths.flatten()
                return [ meta, paths ]
        }

        ch_blastn_combined = ch_combined_blast.join(ch_filtered_combined, remainder: true).map{
            meta, blastn, filteredblastn ->
                if (blastn[0] == null){
                    blastn[0] = []
                }
                if (blastn[1] == null){
                    blastn[1] = []
                }
                if (filteredblastn[0] == null){
                    filteredblastn[0] = []
                }
                if (filteredblastn[1] == null){
                    filteredblastn[1] = []
                }
                return [ meta, blastn[0], blastn[1], filteredblastn[0], filteredblastn[1] ]
            }

        ch_blastn_summary = SUMMARY_BLASTN (
            ch_blastn_combined
        )
        ch_versions = ch_versions.mix(ch_blastn_summary.versions.first())

    // Drop meta of blastn_summary as it is not needed for the combination step of summarizer
        ch_blastn_summary = ch_blastn_summary.summary.map {
                meta, path -> [path]
            }
        }

    //
    // MODULE: Filter out the classified or validated reads
    //
    if (
            (
                ( params.skip_blastn && params.enable_filter ) || params.filter_with_kraken2
            ) && !params.filter_trimmed
        ) {
        ch_kraken2filter = RENAME_FASTQ_HEADERS_PRE.out.fastq
            .join(ISOLATE_IDS_FROM_KRAKEN2_TO_BLASTN.out.classified_ids, by:[0])
        FILTER(
            ch_kraken2filter
        )
        ch_versions = ch_versions.mix(FILTER.out.versions.first())

    } else if (
        params.enable_filter && !params.filter_trimmed
        ) {
        ch_blastn2filter = FILTER_BLASTN_IDENTCOV.out.classified_ids.map {
            meta, path ->
                return [ meta + [ id: meta.id.replaceAll("(_R1|_R2)", "") ], path ]
        }
        .map{
            meta, path -> tuple(groupKey(meta, meta.amount_of_files), path)
        }
        .groupTuple(by:[0])
        ch_combined_short_long_id = RENAME_FASTQ_HEADERS_PRE.out.fastq.map {
            meta, path ->
                return [ meta + [ id: meta.id.replaceAll("(_R1|_R2)", "") ], path ]
        }
        ch_blastnfilter = ch_combined_short_long_id.join(
            ch_blastn2filter, by:[0]
        )
        FILTER(
            ch_blastnfilter
        )
        ch_versions = ch_versions.mix(FILTER.out.versions.first())
    } else if (
        (
            ( params.skip_blastn && params.enable_filter ) || params.filter_with_kraken2
        ) && params.filter_trimmed
    ){
        ch_kraken2filter = FASTP.out.reads
            .join(ISOLATE_IDS_FROM_KRAKEN2_TO_BLASTN.out.classified_ids, by:[0])
        FILTER(
            ch_kraken2filter
        )
        ch_versions = ch_versions.mix(FILTER.out.versions.first())
    } else if (
        params.enable_filter && params.filter_trimmed
    ){
        ch_blastn2filter = FILTER_BLASTN_IDENTCOV.out.classified_ids.map {
            meta, path ->
                return [ meta + [ id: meta.id.replaceAll("(_R1|_R2)", "") ], path ]
        }
        .map{
            meta, path -> tuple(groupKey(meta, meta.amount_of_files), path)
        }
        .groupTuple(by:[0])

        ch_combined_short_long_id = FASTP.out.reads.map {
            meta, path ->
                return [ meta + [ id: meta.id.replaceAll("(_R1|_R2)", "") ], path ]
        }
        ch_blastnfilter = ch_combined_short_long_id.join(
            ch_blastn2filter, by:[0]
        )
        FILTER(
            ch_blastnfilter
        )
        ch_versions = ch_versions.mix(FILTER.out.versions.first())
    }

    //
    // MODULE: Rename headers after filtering
    //
    if ( params.enable_filter ) {
    ch_headers = RENAME_FASTQ_HEADERS_PRE.out.headers.map {
        meta, path ->
            return [ meta + [ id: meta.id.replaceAll("(_R1|_R2)", "") ], path ]
    }

    ch_filtered2rename = FILTER.out.filtered.map {
        meta, path ->
            return [ meta + [ id: meta.id.replaceAll("(_R1|_R2)", "") ], path ]
    }

    ch_rename_filtered = ch_filtered2rename.join(ch_headers, by:[0])

    RENAME_FASTQ_HEADERS_AFTER(
        ch_rename_filtered
    )
    }
    //
    // MODULE: Summarize the classification process
    //
    if (!params.skip_blastn){
    ch_summary = ch_kraken2_summary.mix(ch_blastn_summary).collect().map {
            item -> [['id': "summary_of_kraken2_and_blastn"], item]
        }
    } else {
        ch_summary = ch_kraken2_summary.collect().map {
            item -> [['id': "summary_of_kraken2"], item]
        }
    }

    ch_summary = SUMMARIZER (
        ch_summary
    )
    ch_versions = ch_versions.mix(ch_summary.versions)

    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
