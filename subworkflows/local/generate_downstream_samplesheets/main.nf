//
// Subworkflow with functionality specific to the nf-core/createtaxdb pipeline
//

workflow SAMPLESHEET_TAXPROFILER {
    take:
    ch_reads

    main:
    format = 'csv' // most common format in nf-core

    // Make your samplesheet channel construct here depending on your downstream
    ch_list_for_samplesheet = ch_reads
        .map {
            meta, reads ->
                def out_path            = file(params.outdir).toString() + '/filter/filtered/'
                def sample              = meta.id
                def run_accession       = meta.id - "_longReads"
                def instrument_platform = !meta.long_reads ? "ILLUMINA" : "OXFORD_NANOPORE"
                def fastq_1             = meta.single_end  ? out_path + reads.getName() : out_path + reads[0].getName()
                def fastq_2             = !meta.single_end ? out_path + reads[1].getName() : ""
                def fasta               = ""
            [ sample: sample, run_accession:run_accession, instrument_platform:instrument_platform, fastq_1:fastq_1, fastq_2:fastq_2, fasta:fasta ]
        }

    channelToSamplesheet(ch_list_for_samplesheet,"${params.outdir}/downstream_samplesheets/taxprofiler", format)

}

workflow SAMPLESHEET_MAG {
    //
    // MAG doesn't take PE-data & SE-data in a single samplesheet
    //
    take:
    ch_reads

    main:
    format = 'csv' // most common format in nf-core

    // Combine the short and long reads belonging to the same sample
    ch_reads
        .map{meta, reads ->
                tuple( groupKey(meta.id - "_longReads", 2), meta, reads)
            }
        .groupTuple(remainder: true)
        .map{key, meta, reads ->
                // Pair each read entry with its metadata, then separate short and long reads.
                def short_files = []
                def long_files  = []
                [meta, reads].transpose().each { m, r ->
                    def files = r instanceof List ? r : [r]
                    if (m.long_reads) { long_files.addAll(files) } else { short_files.addAll(files) }
                }
                def new_meta = [
                    id: key,
                    run: key,
                    single_end: meta[0].single_end,
                    long_reads: !long_files.isEmpty(),
                    short_reads: !short_files.isEmpty(),
                    short_reads_paired: short_files.size() > 1
                    ]
                // Short reads precede long reads so the emitter can index them directly.
                [new_meta, short_files + long_files]
            }
        .tap{ ch_reads_grouped }

    // Make your samplesheet channel construct here depending on your downstream
    ch_list_for_samplesheet = ch_reads_grouped
        .map {
            meta, reads ->
                def out_path           = file(params.outdir).toString() + '/filter/filtered/'
                def sample             = meta.id
                def run                = meta.run
                def group              = 0                                                                                     // only used for co-abundance in binning
                def has_short_reads    = meta.short_reads
                def has_long_reads     = meta.long_reads
                def short_reads_1      = has_short_reads ? out_path + reads[0].getName() : ""
                def short_reads_2      = meta.short_reads_paired ? out_path + reads[1].getName() : ""
                def long_reads         = has_long_reads ? out_path + reads.last().getName() : ""                                // If long reads, take final element
                def short_reads_platform = has_short_reads ? "ILLUMINA" : ""
                def long_reads_platform  = has_long_reads ? "OXFORD_NANOPORE" : ""
            [sample: sample, run: run, group: group, short_reads_1: short_reads_1, short_reads_2: short_reads_2, long_reads: long_reads, short_reads_platform: short_reads_platform, long_reads_platform: long_reads_platform]
        }
        .branch{
            se: it.short_reads_2 ==""
            pe: true
        }

    channelToSamplesheet(ch_list_for_samplesheet.pe,"${params.outdir}/downstream_samplesheets/mag-pe", format)
    channelToSamplesheet(ch_list_for_samplesheet.se, "${params.outdir}/downstream_samplesheets/mag-se", format)

}

workflow GENERATE_DOWNSTREAM_SAMPLESHEETS {
    take:
    ch_reads

    main:
    def downstreampipeline_names = params.generate_pipeline_samplesheets.split(",")

    if ( downstreampipeline_names.contains('taxprofiler')) {
        SAMPLESHEET_TAXPROFILER(ch_reads)
    }
    if ( downstreampipeline_names.contains('mag')) {
        SAMPLESHEET_MAG(ch_reads)
    }

}


// Constructs the header string and then the strings of each row, and
def channelToSamplesheet(ch_list_for_samplesheet, path, format) {
    def format_sep = ["csv":",", "tsv":"\t", "txt":"\t"][format]

    def ch_header = ch_list_for_samplesheet

    ch_header
        .first()
        .map{ it.keySet().join(format_sep) }
        .concat( ch_list_for_samplesheet.map{ it.values().join(format_sep) })
        .collectFile(
            name:"${path}.${format}",
            newLine: true,
            sort: false
        )
}
