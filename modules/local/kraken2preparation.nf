process KRAKEN2PREPARATION {
    tag "$meta.id"
    label 'process_high'

    conda "conda-forge::sed=4.8 conda-forge::tar=1.34"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'nf-core/ubuntu:22.04' }"

    input:
    tuple val(meta), path(db)

    output:
    path( "database/" ) , emit: db
    tuple val("${task.process}"), val('tar'), eval('tar --version | grep -oP "tar \\(GNU tar\\) \\K\\d+(\\.\\d+)*"'), emit: versions_tar, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir db_tmp
    tar -xf "${db}" -C db_tmp
    mkdir database
    mv `find db_tmp/ -name "*.k2d"` database/
    """
}
