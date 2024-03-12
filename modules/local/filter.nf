process FILTER {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::seqkit=2.8.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.8.0--h9ee0642_0':
        'biocontainers/seqkit:2.8.0--h9ee0642_0'}"

    input:
    tuple val(meta), path(fastq), path(ids_to_remove)

    output:
    tuple val(meta), path('*.fastq.gz') , emit: filtered
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Extract the sequences from the fastq.gz files
    if [[ "${fastq}" == *" "* ]]; then
        IFS=' ' read -ra array <<< "${fastq}"
        IFS=' ' read -ra array2 <<< "${ids_to_remove}"
            COUNTER=0
        for element in "\${array[@]}"
        do
            COUNTER=\$((COUNTER+1))
            seqkit grep -v -f \${array2[\$(COUNTER-1)]} \$element -o \$(echo ${meta.id})_R\$(echo \$COUNTER)_filtered_renamed.fastq.gz
        done
    else
        seqkit grep -v -f ${ids_to_remove} ${fastq} -o ${meta.id}_filtered_renamed.fastq.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed -E 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/')
    END_VERSIONS
    """
}
