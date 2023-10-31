process KRAKEN2PREPARATION {
    
    input:
    path db 

    output:
    path( "database/" ), emit: db
    path "versions.yml", emit: versions

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