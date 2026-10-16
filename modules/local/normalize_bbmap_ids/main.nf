process NORMALIZE_BBMAP_IDS {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::gawk=5.3.0"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gawk:5.3.0'
        : 'biocontainers/gawk:5.3.0'}"

    input:
    tuple val(meta), path(reads), path(ids)

    output:
    tuple val(meta), path(reads), path('*normalized_ids.txt'), emit: reads_and_ids
    tuple val("${task.process}"), val('gawk'), eval('awk -Wversion | sed "1!d; s/.*Awk //; s/,.*//"'), emit: versions_gawk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # bbmap's filterbyname.sh pairs the input reads and matches the supplied names
    # only against the first mate (R1). Identifiers that originate from the second
    # mate carry a trailing '/2' and therefore never match, so pairs flagged in R2
    # only would not be removed. Rewriting '/(1|2)' to '/1' (and de-duplicating)
    # makes every flagged mate remove the pair from both mates.
    cat ${ids} \\
        | sed -E 's#/[12]\$#/1#' \\
        | awk '!seen[\$0]++' \\
        > ${meta.id}.normalized_ids.txt
    """
}
