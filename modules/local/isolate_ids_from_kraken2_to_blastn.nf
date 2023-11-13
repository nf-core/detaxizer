process ISOLATE_IDS_FROM_KRAKEN2_TO_BLASTN {

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(kraken2results), path(tax2filter)


    output:
    tuple val(meta), path('*classified.txt'), emit: classified
    path "versions.yml", emit: versions

    script:
    """
    while IFS= read -r line
    do
    line=\$(echo -n "\$line" | tr -d '\n')
    if [ \$(grep -c "\$line" $kraken2results) -gt 0 ]; then
        grep "\$line" $kraken2results >> ${meta.id}.classified.txt
    else
        if [ ! -f "${meta.id}.classified.txt" ]; then
            touch ${meta.id}.classified.txt
        fi
    fi
    done < $tax2filter

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    grep: \$(grep --version  | sed -n 1p  | sed 's/grep (GNU grep) //')
    """
}
