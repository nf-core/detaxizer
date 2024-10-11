//
// Subworkflow with functionality specific to the nf-core/createtaxdb pipeline
//

workflow SAMPLESHEET_TAXPROFILER{
    take:
    ch_reads

    main:
    ch_header  = Channel.empty()
    format     = 'csv' // most common format in nf-core
    format_sep = ','

    // Make your samplesheet channel construct here depending on your downstream
    ch_list_for_samplesheet = ch_reads
                                .map {
                                    meta, reads ->
                                        def out_path            = file(params.outdir).toString() + '/filter/filtered/'
                                        def sample              = meta.id
                                        def run_accession       = meta.id - "_longReads"
                                        def instrument_platform = !meta.long_reads ? "ILLUMINA" : "OXFORD_NANOPORE"
                                        def fastq_1             = meta.single_end  ?  out_path + reads.getName(): out_path + reads[0].getName()
                                        def fastq_2             = !meta.single_end ? out_path + reads[1].getName() : ""
                                        def fasta               = ""
                                    [ sample: sample, run_accession:run_accession, instrument_platform:instrument_platform, fastq_1:fastq_1, fastq_2:fastq_2, fasta:fasta ]
                                }
                                .tap{ ch_colnames } //ch_header exists

    channelToSamplesheet(ch_colnames, ch_list_for_samplesheet, format, format_sep)

}

workflow GENERATE_DOWNSTREAM_SAMPLESHEETS {
    take:
    ch_reads

    main:
    def downstreampipeline_names = params.generate_pipeline_samplesheets.split(",")

    if ( downstreampipeline_names.contains('taxprofiler')) {
        SAMPLESHEET_TAXPROFILER(ch_reads)
    } else {
        error "Unsupported downstream pipeline: ${params.generate_pipeline_samplesheets}"
    }
}

def channelToSamplesheet(ch_header,ch_list_for_samplesheet, format, format_sep) {
    // Constructs the header string and then the strings of each row, and
    // finally concatenates for saving.
    ch_header
        .first()
        .map{ it.keySet().join(format_sep) }
        .concat( ch_list_for_samplesheet.map{ it.values().join(format_sep) })
        .collectFile(
            name:"${params.outdir}/downstream_samplesheet/${params.generate_pipeline_samplesheets}.${format}",
            newLine: true,
            sort: false
        )
}
