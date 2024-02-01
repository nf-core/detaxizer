process PREPARE_FASTA4BLASTN {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::seqkit=2.6.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit%3A2.6.0--h9ee0642_0':
        'biocontainers/seqkit:2.6.0--h9ee0642_0'}"

    input:
    tuple val(meta), path(trimmedreads), path(kraken2results)

    output:
    tuple val(meta), path("*.fa.gz"), emit: fasta
    path("versions.yml")            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    if [ "$meta.single_end" == "true" ]; then
        seqkit fq2fa ${trimmedreads} -o out.fa.gz
        seqkit grep -f ${kraken2results} out.fa.gz -o ${meta.id}.fa.gz
        rm out.fa.gz
    else
        seqkit fq2fa ${trimmedreads[0]} -o out.fa.gz
        seqkit grep -f ${kraken2results} out.fa.gz -o ${meta.id}_R1.fa.gz
        rm out.fa.gz
        seqkit fq2fa ${trimmedreads[1]} -o out.fa.gz
        seqkit grep -f ${kraken2results} out.fa.gz -o ${meta.id}_R2.fa.gz
        rm out.fa.gz
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed -E 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/')
    END_VERSIONS
    """
}
