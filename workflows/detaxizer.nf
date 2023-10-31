/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowDetaxizer.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl=2

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { KRAKEN2PREPARATION } from '../modules/local/kraken2preparation'
include { PARSE_KRAKEN2REPORT } from '../modules/local/parse_kraken2report'
include { ISOLATE_IDS_FROM_KRAKEN2_TO_BLASTN } from '../modules/local/isolate_ids_from_kraken2_to_blastn'
include { PREPARE_FASTA4BLASTN } from '../modules/local/prepare_fasta4blastn'
include { FILTER_BLASTN_IDENTCOV } from '../modules/local/filter_blastn_identcov'
include { SUMMARY_KRAKEN2 } from '../modules/local/summary_kraken2'
include { SUMMARY_BLASTN } from '../modules/local/summary_blastn'
include { SUMMARIZER } from '../modules/local/summarizer'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { FASTP } from '../modules/nf-core/fastp/main'
include { KRAKEN2_KRAKEN2 } from '../modules/nf-core/kraken2/kraken2/main'
include { BLAST_BLASTN } from '../modules/nf-core/blast/blastn/main'
include { BLAST_MAKEBLASTDB } from '../modules/nf-core/blast/makeblastdb'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []




workflow DETAXIZER {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    //INPUT_CHECK (
    //    file(params.input)
    //)
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema
    // TODO: Make the map for single ended 
    ch_input = Channel.fromSamplesheet('input')


    // check whether the sample sheet is correctly formated    
    ch_input.map {
        meta, fastq_1, fastq_2, fastq_3 ->
        if (!fastq_1 && !fastq_3){
            error("Please provide at least one single end file as input in the sample sheet for ${meta.id}.")
        } else if (!fastq_1 && fastq_2 && fastq_3){
            error("Please provide single end reads in following format in the sample sheet: base name, fastq_1,,fastq_3. fastq_1 is the short read file, fastq_3 the long read file. The wrongly formated entry is ${meta.id}.")
        }
    }

  

    ch_input.branch { 
	    shortReads: it[1]
	    }.set { 
            ch_short 
        }

    ch_short.shortReads.map{
        meta, fastq_1, fastq_2, fastq_3 ->
        if (fastq_2){
            def newMeta = meta.clone()
            newMeta.single_end = false
            newMeta.long_reads = false
            return [newMeta, [fastq_1, fastq_2]]
        } else {
            def newMeta = meta.clone()
            newMeta.id = "${newMeta.id}_R1"
            newMeta.single_end = true
            newMeta.long_reads = false
            return [newMeta, fastq_1]
        }
    }.set{
        ch_short
    }

    ch_input.branch {
        longReads: it[3]
    }.set { 
        ch_long
    }

    ch_long.longReads.map {
        meta, fastq_1, fastq_2, fastq_3 ->
        def newMeta = meta.clone()
        newMeta.id = "${newMeta.id}_R3"
        newMeta.single_end = true
        newMeta.long_reads = true
        return [newMeta, fastq_3]
    }.set {
        ch_long
    }


    ch_combined_short_long = ch_short.mix(ch_long)


    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_combined_short_long
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Run fastp
    //

    FASTP (
        ch_combined_short_long,
        [],
        params.fastp_save_trimmed_fail,
        []
    )
    ch_versions = ch_versions.mix(FASTP.out.versions)

    //
    // MODULE: Prepare Kraken2 Database
    //

    KRAKEN2PREPARATION (
        params.kraken2db
    )
    ch_versions = ch_versions.mix(KRAKEN2PREPARATION.out.versions)

    // 
    // MODULE: Run Kraken2
    //

    KRAKEN2_KRAKEN2 (
        FASTP.out.reads,
        KRAKEN2PREPARATION.out.db,
        params.save_output_fastqs,
        params.save_reads_assignment
    )
    ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions)

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
    
    KRAKEN2_KRAKEN2.out.classified_reads_assignment.combine(PARSE_KRAKEN2REPORT.out.txt).set{ ch_combined }

    ISOLATE_IDS_FROM_KRAKEN2_TO_BLASTN (
        ch_combined
    )

    ch_versions = ch_versions.mix(ISOLATE_IDS_FROM_KRAKEN2_TO_BLASTN.out.versions)

    //
    // MODULE: Summarize the kraken2 results and the isolated kraken2 hits
    //

    // preparation of the channels
    ch_prepare_summary_kraken2 = KRAKEN2_KRAKEN2.out.classified_reads_assignment.join(ISOLATE_IDS_FROM_KRAKEN2_TO_BLASTN.out.classified)
    ch_prepare_summary_kraken2 = ch_prepare_summary_kraken2.map { 
        meta, path1, path2 -> 
        [ ['id': meta.id ], [ path1, path2 ] ] 
    }

    ch_combined_kraken2 = ch_prepare_summary_kraken2
        .map { meta, path -> 
        [ ['id': meta.id.replaceAll("(_R1|_R2|_R3)", "")], path ]
        }
        .groupTuple(by: [0])
        .map { meta, path ->
            path = path.flatten()
            return [meta, path]
            }
        .map { meta, path ->
            if(path.size() == 4) {
                def newMeta = meta.clone()
                newMeta.short_and_long_reads = true
                return [newMeta, path]
            } else {
                def newMeta = meta.clone()
                newMeta.short_and_long_reads = false
                return [newMeta, path]
            }
        }

    // run of the process
    ch_kraken2_summary = SUMMARY_KRAKEN2(
        ch_combined_kraken2
        ).map {
            meta, path -> [path]
        }
    //
    // MODULE: Extract the hits to fasta format
    //
    // TODO: skip from here onward if isolate from kraken2 is empyty using .branch
    ch_combined = FASTP.out.reads
        .join(
            ISOLATE_IDS_FROM_KRAKEN2_TO_BLASTN.out.classified, by: [0]
        ) 
    // TODO: REPLACE awk AND seqtk WITH TOOLS THAT CAN DISPLAY THEIR VERSION NUMBER
    PREPARE_FASTA4BLASTN (
        ch_combined
    )
    ch_versions = ch_versions.mix(PREPARE_FASTA4BLASTN.out.versions)    

    // 
    // MODULE: Run BLASTN if --skip_blastn = false
    //
    
    ch_reference_fasta = Channel.empty()
    
    if (!params.skip_blastn) {  // If skip_blastn is false, then execute the process
        ch_reference_fasta =  file( params.fasta )
    }
    
    BLAST_MAKEBLASTDB (
            ch_reference_fasta
    )
    ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)

    //PREPARE_FASTA4BLASTN.out.fasta.dump(tag: "PREPARE_FASTA4BLASTN")
    ch_fasta4blastn = PREPARE_FASTA4BLASTN.out.fasta
        .flatMap { meta, fastaList ->
            if (fastaList.size() == 2) {
        // Handle paired-end FASTA files
            return [
                [['id': "${meta.id}_R1", 'single_end': false], fastaList[0]],
                [['id': "${meta.id}_R2", 'single_end': false], fastaList[1]]
            ]

        } else {
        // Handle single-end FASTA files
        return [ [ [ 'id': "${meta.id}", 'single_end': true], fastaList ] ]
        } 
      
        }
    //ch_fasta4blastn.dump(tag: "ch_fasta4blastn")
    
    BLAST_BLASTN (
        ch_fasta4blastn,
        BLAST_MAKEBLASTDB.out.db
    )
    
    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions)
    
    ch_combined_blast = BLAST_BLASTN.out.txt
                                .map { meta, path -> 
        [ ['id': meta.id.replaceAll("(_R1|_R2|_R3)", "")], path ]
    }


    ch_combined_blast = ch_combined_blast
                               .groupTuple(by: [0])
                               .map { meta, paths -> [ meta, paths.flatten() ] }

    FILTER_BLASTN_IDENTCOV (
        BLAST_BLASTN.out.txt
    )
    ch_versions = ch_versions.mix(FILTER_BLASTN_IDENTCOV.out.versions)
    
    ch_filtered_combined = FILTER_BLASTN_IDENTCOV.out.classified.map { 
        meta, path -> 
        [ ['id': meta.id.replaceAll("(_R1|_R2|_R3)", "")], path ]
    }
    ch_filtered_combined = ch_filtered_combined
        .groupTuple ()
        .map {
            meta, paths -> 
            paths = paths.flatten()
            return [meta, paths] 
        }
    
    ch_blastn_combined = ch_combined_blast.join(ch_filtered_combined, remainder: true)
    //ch_blastn_combined.dump()
    ch_blastn_combined = ch_blastn_combined.map { 
        meta, blastn, filteredblastn -> 
            if (blastn[0] == null){
                blastn[0] = "${projectDir}/assets/NO_FILE1"
            }
            if (blastn[1] == null){
                blastn[1] = "${projectDir}/assets/NO_FILE2"
            }
            if (blastn[2] == null){
                blastn[2] = "${projectDir}/assets/NO_FILE3"
            }
            if (filteredblastn[0] == null){
                filteredblastn[0] = "${projectDir}/assets/NO_FILE4"
            }
            if (filteredblastn[1] == null){
                filteredblastn[1] = "${projectDir}/assets/NO_FILE5"
            }
            if (filteredblastn[2] == null){
                filteredblastn[2] = "${projectDir}/assets/NO_FILE6"
            }
            return [meta, blastn[0], blastn[1], blastn[2], filteredblastn[0], filteredblastn[1], filteredblastn[2]]
        }

    // TODO: Get Software version and prepare it for multiqc
    ch_blastn_summary = SUMMARY_BLASTN (
        ch_blastn_combined
    )

    //
    // MODULE: Summarize the classification process
    //

    //ch_kraken2_summary = ch_kraken2_summary.map { meta, paths -> [paths] }
    ch_blastn_summary = ch_blastn_summary.map { meta, paths -> [paths]}

    ch_summary = ch_kraken2_summary.mix(ch_blastn_summary).collect()

    SUMMARIZER (
        ch_summary
    )

    

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowDetaxizer.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowDetaxizer.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()

    
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
