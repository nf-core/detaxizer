process RENAME_FASTQ_HEADERS_PRE {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.12.0 numpy=1.26.3 dnaio=1.2.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dnaio:1.2.2--py312hf67a6ed_0' :
        'quay.io/biocontainers/dnaio:1.2.2--py312hf67a6ed_0' }"

    input:
    tuple val(meta), path(inputfastq)

    output:
    tuple val(meta), path('*_headers*.txt.gz')  , emit: headers
    tuple val(meta), path('*.fastq.gz')         , emit: fastq
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    rename_fastq_headers_pre.py -i $inputfastq -o $meta.id

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
    END_VERSIONS
    """
}
