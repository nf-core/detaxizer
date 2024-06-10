process PARSE_KRAKEN2REPORT {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(kraken2report)

    output:
    tuple val(meta), path ("taxa_to_filter.txt"), emit: to_filter
    path "versions.yml",                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    parse_kraken2report.py -i $kraken2report -t $params.tax2filter

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
