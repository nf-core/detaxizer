process PARSE_KRAKEN2REPORT {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.12.2"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(kraken2report)

    output:
    tuple val(meta), path ("taxa_to_filter.txt"), emit: to_filter
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //g"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    parse_kraken2report.py -i $kraken2report -t "$params.tax2filter"
    """
}
