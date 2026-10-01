process RENAME_FASTQ_HEADERS_PRE {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.12.6 numpy=1.26.3 dnaio=1.2.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dnaio:1.2.2--py312hf67a6ed_0' :
        'quay.io/biocontainers/dnaio:1.2.2--py312hf67a6ed_0' }"

    input:
    tuple val(meta), path(inputfastq)

    output:
    tuple val(meta), path('*_headers*.txt.gz')  , emit: headers
    tuple val(meta), path('*.fastq.gz')         , emit: fastq
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //g"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    rename_fastq_headers_pre.py -i $inputfastq -o $meta.id
    """
}
