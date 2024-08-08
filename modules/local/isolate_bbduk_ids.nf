process ISOLATE_BBDUK_IDS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::seqkit=2.8.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.8.0--h9ee0642_0':
        'biocontainers/seqkit:2.8.0--h9ee0642_0'}"

    input:
    tuple val(meta), path(contamination)

    output:
    tuple val(meta), path('*.bbduk.txt')    , emit: classified_ids
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    if [ "$meta.single_end" == "true" ]; then
        seqkit seq -n $contamination > ${meta.id}.bbduk.txt
    else
        seqkit seq -n ${contamination[0]} > read_ids1.txt
        seqkit seq -n ${contamination[1]} > read_ids2.txt
        awk '!seen[\$0]++' read_ids1.txt read_ids2.txt > ${meta.id}.bbduk.txt
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed -E 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/')
    END_VERSIONS
    """
}
