/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                                                    } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                                   } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                                          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                                      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                                    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                                    } from '../subworkflows/local/utils_nfcore_detaxizer_pipeline'
include { getGenomeAttribute                                        } from '../subworkflows/local/utils_nfcore_detaxizer_pipeline'
include { GENERATE_DOWNSTREAM_SAMPLESHEETS                          } from '../subworkflows/local/generate_downstream_samplesheets/main.nf'

include { FASTP                                                     } from '../modules/nf-core/fastp/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_KRAKEN2                        } from '../modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_POST_CLASSIFICATION_FILTERED   } from '../modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_POST_CLASSIFICATION_REMOVED    } from '../modules/nf-core/kraken2/kraken2/main'
include { BBMAP_BBDUK                                               } from '../modules/nf-core/bbmap/bbduk/main'
include { BLAST_BLASTN                                              } from '../modules/nf-core/blast/blastn/main'
include { BLAST_MAKEBLASTDB                                         } from '../modules/nf-core/blast/makeblastdb/main'
include { BBMAP_FILTERBYNAME                                        } from '../modules/nf-core/bbmap/filterbyname/main'
include { BBMAP_FILTERBYNAME as BBMAP_FILTERBYNAME_REMOVED          } from '../modules/nf-core/bbmap/filterbyname/main'

include { RENAME_FASTQ_HEADERS_PRE                                  } from '../modules/local/rename_fastq_headers_pre'
include { KRAKEN2PREPARATION                                        } from '../modules/local/kraken2preparation'
include { PARSE_KRAKEN2REPORT                                       } from '../modules/local/parse_kraken2report'
include { ISOLATE_KRAKEN2_IDS                                       } from '../modules/local/isolate_kraken2_ids'
include { MAP_KRAKEN2SEQIDS_TO_FQHEADERS                            } from '../modules/local/map_kraken2seqids_to_fqheaders'
include { ISOLATE_BBDUK_IDS                                         } from '../modules/local/isolate_bbduk_ids'
include { MERGE_IDS                                                 } from '../modules/local/merge_ids'
include { PREPARE_FASTA4BLASTN                                      } from '../modules/local/prepare_fasta4blastn'
include { FILTER_BLASTN_IDENTCOV                                    } from '../modules/local/filter_blastn_identcov'
include { FILTER                                                    } from '../modules/local/filter'
include { RENAME_FASTQ_HEADERS_AFTER                                } from '../modules/local/rename_fastq_headers_after'
include { SUMMARY_CLASSIFICATION                                    } from '../modules/local/summary_classification'
include { SUMMARY_BLASTN                                            } from '../modules/local/summary_blastn'
include { SUMMARIZER                                                } from '../modules/local/summarizer'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// specify the ch_fasta_blastn channel if it is not provided via --fasta_blastn
def ch_fasta_blastn = Channel.empty()

if ( !params.fasta_blastn && params.validation_blastn ) {
    ch_fasta_blastn = Channel.fromPath(getGenomeAttribute('fasta'))
} else if ( params.validation_blastn ){
    // If params.fasta_blastn is there, use it for the creation of the blastn database
    ch_fasta_blastn = Channel.fromPath(params.fasta_blastn)
}

// specify the ch_fasta_bbduk channel if it is not provided via --fasta_bbduk

def ch_fasta_bbduk = Channel.empty()

if ( !params.fasta_bbduk && params.classification_bbduk ) {
    ch_fasta_bbduk = Channel.fromPath(getGenomeAttribute('fasta'))
} else if ( params.classification_bbduk ){
    // If params.fasta_bbduk is there, use it for the creation of the blastn database
    ch_fasta_bbduk = Channel.fromPath(params.fasta_bbduk)
}

workflow NFCORE_DETAXIZER {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_filtered_reads = Channel.empty()
    ch_removed_reads  = Channel.empty()

    ch_short = ch_samplesheet.branch {
        shortReads: it[1]
        }.shortReads.map{
        meta, short_reads_fastq_1, short_reads_fastq_2, long_reads_fastq_1 ->
            if (short_reads_fastq_2){
                return [meta + [ single_end: false, long_reads: false , amount_of_files: 2 ], [ short_reads_fastq_1, short_reads_fastq_2 ] ]
            } else {
                return [meta + [ id: "${meta.id}_R1", single_end: true, long_reads: false, amount_of_files: 1 ], short_reads_fastq_1 ]
            }
    }

    ch_long = ch_samplesheet.branch {
        longReads: it[3]
    }.longReads.map {
        meta, short_reads_fastq_1, short_reads_fastq_2, long_reads_fastq_1 ->
            return [meta + [ id: "${meta.id}_longReads", single_end: true, long_reads: true, amount_of_files: 1 ], long_reads_fastq_1 ]
    }

    ch_short_long = ch_short.mix(ch_long)

    //
    // MODULE: Rename Fastq headers
    //
    if ( params.filtering_tool == 'seqkit' ) {
        RENAME_FASTQ_HEADERS_PRE(ch_short_long)
        ch_fastq_input = RENAME_FASTQ_HEADERS_PRE.out.fastq
    } else {
        ch_fastq_input = ch_short_long
    }

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_fastq_input
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Run fastp
    //
    if (params.filter_trimmed || params.preprocessing) {

    FASTP (
        ch_fastq_input,
        [],
        [],
        params.fastp_save_trimmed_fail,
        []
    )

    ch_fastq_for_classification = FASTP.out.reads
    ch_versions = ch_versions.mix(FASTP.out.versions.first())
    } else {
        ch_fastq_for_classification = ch_fastq_input
    }
    //////////////////////////////////////////////////
    //  Classification
    //////////////////////////////////////////////////

    if ( params.classification_kraken2_post_filtering || (!params.classification_kraken2 && !params.classification_bbduk) || (params.classification_kraken2) ){
        //
        // MODULE: Prepare Kraken2 Database
        //
        ch_kraken2_db = Channel.fromPath(params.kraken2db).map {
                item -> [['id': "kraken2_db"], item]
            }

        KRAKEN2PREPARATION (
            ch_kraken2_db
        )
        ch_versions = ch_versions.mix(KRAKEN2PREPARATION.out.versions.first())
    }


    if ((!params.classification_bbduk && !params.classification_kraken2) || (params.classification_kraken2)) {

        //
        // MODULE: Run Kraken2
        //
        KRAKEN2_KRAKEN2 (
            ch_fastq_for_classification,
            KRAKEN2PREPARATION.out.db.first(),
            params.save_output_fastqs,
            true
        )
        ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions.first())

        if ( params.filtering_tool == 'bbmap' ) {
            MAP_KRAKEN2SEQIDS_TO_FQHEADERS(
                ch_fastq_for_classification
                    .join(KRAKEN2_KRAKEN2.out.classified_reads_assignment, by: 0)
                    .map { meta, reads, classification -> [meta, reads, classification] }
            )
            ch_versions = ch_versions.mix(MAP_KRAKEN2SEQIDS_TO_FQHEADERS.out.versions.first())
        }

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

        if ( params.filtering_tool == 'bbmap' ) {
            ch_combined = KRAKEN2_KRAKEN2.out.classified_reads_assignment
                .combine(ch_parsed_kraken2_report)
                .join(MAP_KRAKEN2SEQIDS_TO_FQHEADERS.out.mapping, by:[0])
                .map { meta, kraken2results, tax, mapping -> [meta, kraken2results, tax, mapping] }
        } else {
            ch_combined = KRAKEN2_KRAKEN2.out.classified_reads_assignment
                .combine(ch_parsed_kraken2_report)
                .map { meta, kraken2results, tax -> [meta, kraken2results, tax, null] }
        }

        ISOLATE_KRAKEN2_IDS (
            ch_combined
        )

        ch_versions = ch_versions.mix(ISOLATE_KRAKEN2_IDS.out.versions.first())

        }

    if (params.classification_bbduk) {

        //
        // MODULE: Run bbduk
        //
        BBMAP_BBDUK (
            ch_fastq_for_classification,
            ch_fasta_bbduk.first()
        )
        ch_versions = ch_versions.mix(BBMAP_BBDUK.out.versions.first())

        //
        // MODULE: Run ISOLATE_BBDUK_IDS
        //
        ISOLATE_BBDUK_IDS(
            BBMAP_BBDUK.out.contaminated_reads
        )
        ch_versions = ch_versions.mix(ISOLATE_BBDUK_IDS.out.versions.first())


    }

    // Prepare MERGE_IDS Channel (with or without merging of IDs)

    if (params.classification_bbduk && params.classification_kraken2){

        //
        // MODULE: Merge IDs
        //
        MERGE_IDS(
            ISOLATE_KRAKEN2_IDS.out.classified_ids.join(
                ISOLATE_BBDUK_IDS.out.classified_ids, by: [0]
            ).map{
                meta, path1, path2 ->
                    [meta,[path1,path2]]
            }
        )


    } else if (params.classification_bbduk && !params.classification_kraken2){

        //
        // MODULE: Merge IDs
        //
        MERGE_IDS(
            ISOLATE_BBDUK_IDS.out.classified_ids
        )


    } else if (params.classification_kraken2 || (!params.classification_kraken2 && !params.classification_bbduk)){

        //
        // MODULE: Merge IDs
        //
        MERGE_IDS(
            ISOLATE_KRAKEN2_IDS.out.classified_ids
        )

    }

    ch_versions = ch_versions.mix(MERGE_IDS.out.versions.first())

    //
    // MODULE: Summarize the classification results
    //

    SUMMARY_CLASSIFICATION(
        MERGE_IDS.out.classified_ids
    )

    // Drop meta of kraken2_summary as it is not needed for the combination step of summarizer
    ch_classification_summary = SUMMARY_CLASSIFICATION.out.summary.map {
            meta, path -> [path]
    }
    ch_versions = ch_versions.mix(SUMMARY_CLASSIFICATION.out.versions.first())

    //////////////////////////////////////////////////
    //  Validation
    //////////////////////////////////////////////////

    if (params.validation_blastn) {

        //
        // MODULE: Extract the hits to fasta format
        //
        ch_combined = ch_fastq_for_classification
        .join(
            MERGE_IDS.out.classified_ids, by: [0]
        )


        PREPARE_FASTA4BLASTN (
            ch_combined
        )

        ch_versions = ch_versions.mix(PREPARE_FASTA4BLASTN.out.versions.first())

        //
        // MODULE: Run BLASTN
        //
        ch_reference_fasta = ch_fasta_blastn

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
    if ( !params.skip_filter ) {

        ch_reads_for_filter = params.filter_trimmed ? ch_fastq_for_classification : ch_fastq_input

        if ( params.validation_blastn && !params.filter_with_classification ) {

            ch_blastn2filter = FILTER_BLASTN_IDENTCOV.out.classified_ids.map { meta, path ->
                [ meta + [ id: meta.id.replaceAll("(_R1|_R2)", "") ], path ]
            }
            .map { meta, path -> tuple(groupKey(meta, meta.amount_of_files), path) }
            .groupTuple(by:[0])

            ch_reads_with_id = ch_reads_for_filter.map { meta, path ->
                [ meta + [ id: meta.id.replaceAll("(_R1|_R2)", "") ], path ]
            }

            ch_to_filter = ch_reads_with_id.join(ch_blastn2filter, by:[0])

        } else {

            ch_to_filter = ch_reads_for_filter.join(MERGE_IDS.out.classified_ids, by:[0])

        }

        if ( params.filtering_tool == 'seqkit' ) {
            FILTER(
                ch_to_filter
            )
            ch_versions = ch_versions.mix(FILTER.out.versions.first())
            ch_filter_filtered = FILTER.out.filtered
            ch_filter_removed  = FILTER.out.removed
        } else {
            BBMAP_FILTERBYNAME(
                ch_to_filter.map { meta, reads, ids -> tuple(meta, reads) },
                ch_to_filter.map { meta, reads, ids -> ids.toString() },
                Channel.value('fastq.gz'),
                Channel.value(false)
            )
            ch_versions = ch_versions.mix(BBMAP_FILTERBYNAME.out.versions.first())
            ch_filter_filtered = BBMAP_FILTERBYNAME.out.reads
            if ( params.output_removed_reads ) {
                BBMAP_FILTERBYNAME_REMOVED(
                    ch_to_filter.map { meta, reads, ids -> tuple(meta, reads) },
                    ch_to_filter.map { meta, reads, ids -> ids.toString() },
                    Channel.value('fastq.gz'),
                    Channel.value(false)
                )
                ch_versions = ch_versions.mix(BBMAP_FILTERBYNAME_REMOVED.out.versions.first())
                ch_filter_removed = BBMAP_FILTERBYNAME_REMOVED.out.reads
            } else {
                ch_filter_removed = Channel.empty()
            }
        }
    }
    //
    // MODULE: Rename headers after filtering
    //
    if ( !params.skip_filter ) {

        if ( params.filtering_tool == 'seqkit' ) {
            ch_headers = RENAME_FASTQ_HEADERS_PRE.out.headers.map { meta, path ->
                [ meta + [ id: meta.id.replaceAll("(_R1|_R2)", "") ], path ]
            }

            ch_filtered2rename = ch_filter_filtered.map { meta, path ->
                [ meta + [ id: meta.id.replaceAll("(_R1|_R2)", "") ], path ]
            }

            ch_removed2rename = Channel.empty()
            if ( params.output_removed_reads ) {
                ch_removed2rename = ch_filter_removed.map { meta, path ->
                    [ meta + [ id: meta.id.replaceAll("(_R1|_R2)", "") ], path ]
                }
            }

            ch_rename_filtered = ch_filtered2rename.join(ch_headers, by:[0])
            ch_removed2rename = ch_removed2rename.ifEmpty(['empty', []])

            if ( params.output_removed_reads ) {
                RENAME_FASTQ_HEADERS_AFTER(
                    ch_rename_filtered,
                    ch_removed2rename
                )
            } else {
                RENAME_FASTQ_HEADERS_AFTER(
                    ch_rename_filtered,
                    ch_removed2rename.first()
                )
            }
            ch_versions = ch_versions.mix(RENAME_FASTQ_HEADERS_AFTER.out.versions.first())
            ch_filtered_reads = RENAME_FASTQ_HEADERS_AFTER.out.fastq
            ch_removed_reads  = params.output_removed_reads ? RENAME_FASTQ_HEADERS_AFTER.out.fastq_removed : Channel.empty()
        } else {
            ch_filtered_reads = ch_filter_filtered.map { meta, path ->
                [ meta + [ id: meta.id.replaceAll("(_R1|_R2)", "") ], path ]
            }
            ch_removed_reads = Channel.empty()
            if ( params.output_removed_reads ) {
                ch_removed_reads = ch_filter_removed.map { meta, path ->
                    [ meta + [ id: meta.id.replaceAll("(_R1|_R2)", "") ], path ]
                }
            }
        }

        if ( params.classification_kraken2_post_filtering ) {

            KRAKEN2_POST_CLASSIFICATION_FILTERED (
                ch_filtered_reads,
                KRAKEN2PREPARATION.out.db.first(),
                params.save_output_fastqs_filtered,
                true
                )

            ch_versions = ch_versions.mix(KRAKEN2_POST_CLASSIFICATION_FILTERED.out.versions.first())

            if (params.output_removed_reads) {

                KRAKEN2_POST_CLASSIFICATION_REMOVED (
                    ch_removed_reads,
                    KRAKEN2PREPARATION.out.db.first(),
                    params.save_output_fastqs_removed,
                    true
                    )

                ch_versions = ch_versions.mix(KRAKEN2_POST_CLASSIFICATION_REMOVED.out.versions.first())

            }

        }
    }

    //
    // MODULE: Summarize the classification process
    //
    if (params.validation_blastn){

    ch_summary = ch_classification_summary.mix(ch_blastn_summary).collect().map {
            item -> [['id': "summary_of_classification_and_blastn"], item]
        }

    } else {

        ch_summary = ch_classification_summary.collect().map {
            item -> [['id': "summary_of_classification"], item]
        }

    }

    ch_summary = SUMMARIZER (
        ch_summary
    )

    ch_versions = ch_versions.mix(ch_summary.versions)

    if ( params.generate_downstream_samplesheets ) {

        GENERATE_DOWNSTREAM_SAMPLESHEETS ( ch_filtered_reads )

    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'detaxizer_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
