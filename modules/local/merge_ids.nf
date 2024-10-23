process MERGE_IDS {
    tag "$meta.id"
    label 'process_high'

    conda "conda-forge::gawk=5.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'biocontainers/gawk:5.3.0' }"

    input:
    tuple val(meta), path(ids)

    output:
    tuple val(meta), path('*ids.txt')       , emit: classified_ids
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    stringarray=($ids)
    if [ "\${#stringarray[@]}" == 1 ]; then
        cat \${stringarray[0]} > ${meta.id}.ids.txt
    else
        awk '!seen[\$0]++' \${stringarray[0]} \${stringarray[1]} > ${meta.id}.ids.txt
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
