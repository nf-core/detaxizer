process PREPARE_FASTA4BLASTN {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::seqkit=2.8.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.8.0--h9ee0642_0':
        'biocontainers/seqkit:2.8.0--h9ee0642_0'}"

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
        seqkit grep -f ${kraken2results} ${trimmedreads} -o out.fq.gz
        seqkit fq2fa out.fq.gz -o ${meta.id}.fa.gz
        rm out.fq.gz
    else
        seqkit grep -f ${kraken2results} ${trimmedreads[0]} -o out.fq.gz
        seqkit fq2fa out.fq.gz -o ${meta.id}_R1.fa.gz
        rm out.fq.gz
        seqkit grep -f ${kraken2results} ${trimmedreads[1]} -o out.fq.gz
        seqkit fq2fa out.fq.gz -o ${meta.id}_R2.fa.gz
        rm out.fq.gz
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed -E 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/')
    END_VERSIONS
    """
}
