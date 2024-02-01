process KRAKEN2PREPARATION {
    tag "$meta.id"
    label 'process_high'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(db)

    output:
    path( "database/" ) , emit: db
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir db_tmp
    tar -xf "${db}" -C db_tmp
    mkdir database
    mv `find db_tmp/ -name "*.k2d"` database/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tar: \$(tar --version | grep -oP 'tar \\(GNU tar\\) \\K\\d+(\\.\\d+)*')
    END_VERSIONS
    """
}
