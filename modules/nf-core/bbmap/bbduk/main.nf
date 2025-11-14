process BBMAP_BBDUK {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5aae5977ff9de3e01ff962dc495bfa23f4304c676446b5fdf2de5c7edfa2dc4e/data' :
        'community.wave.seqera.io/library/bbmap_pigz:07416fe99b090fa9' }"

    input:
    tuple val(meta), path(reads)
    path contaminants

    output:
    tuple val(meta), path('*.uncontaminated.fastq.gz')  , emit: reads
    tuple val(meta), path('*.contaminated.fastq.gz')    , emit: contaminated_reads
    tuple val(meta), path('*.log')                      , emit: log
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def raw      = meta.single_end ? "in=${reads[0]}" : "in1=${reads[0]} in2=${reads[1]}"
    def trimmed  = meta.single_end ? "out=${prefix}.uncontaminated.fastq.gz" : "out1=${prefix}_1.uncontaminated.fastq.gz out2=${prefix}_2.uncontaminated.fastq.gz"
    def contaminated_reads = meta.single_end ? "outm=${prefix}.contaminated.fastq.gz" : "outm=${prefix}_1.contaminated.fastq.gz outm2=${prefix}_2.contaminated.fastq.gz"
    def contaminants_fa = contaminants ? "ref=$contaminants" : ''
    """
    bbduk.sh \\
        -Xmx${task.memory.toGiga()}g \\
        $raw \\
        $trimmed \\
        $contaminated_reads \\
        threads=$task.cpus \\
        $args \\
        $contaminants_fa \\
        &> ${prefix}.bbduk.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "]")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_command  = meta.single_end ? "echo '' | gzip > ${prefix}.uncontaminated.fastq.gz ; echo '' | gzip > ${prefix}.contaminated.fastq.gz" : "echo '' | gzip > ${prefix}_1.uncontaminated.fastq.gz ; echo '' | gzip > ${prefix}_2.uncontaminated.fastq.gz ; echo '' | gzip > ${prefix}_1.contaminated.fastq.gz ; echo '' | gzip > ${prefix}_2.contaminated.fastq.gz"
    """
    touch ${prefix}.bbduk.log
    $output_command

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "]")
    END_VERSIONS
    """
}
