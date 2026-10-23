process NORMALIZE_BBMAP_IDS {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::gawk=5.3.0"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gawk:5.3.0'
        : 'biocontainers/gawk:5.3.0'}"

    input:
    tuple val(meta), path(ids)

    output:
    tuple val(meta), path('*normalized_ids.txt'), emit: ids
    tuple val("${task.process}"), val('gawk'), eval('awk -Wversion | sed "1!d; s/.*Awk //; s/,.*//"'), emit: versions_gawk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def id_files = [ids].flatten().collect { id_file -> "\"${id_file}\"" }.join(' ')
    """
    normalize_bbmap_ids.sh -o "${meta.id}.normalized_ids.txt" ${id_files}
    """
}
