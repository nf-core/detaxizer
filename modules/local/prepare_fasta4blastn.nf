process PREPARE_FASTA4BLASTN {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::seqkit=2.6.1-0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit%3A2.6.0--h9ee0642_0':
        'biocontainers/seqkit:2.6.0--h9ee0642_0'}"

    input:
    tuple val(meta), path(trimmedreads), path(kraken2results)

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: 2.6.0
    END_VERSIONS

    if [ "$meta.single_end" == "true" ]; then
        seqkit fq2fa ${trimmedreads} -o out.fasta
        seqkit grep -f ${kraken2results} out.fasta -o ${meta.id}.fasta
        rm out.fasta
    else
        seqkit fq2fa ${trimmedreads[0]} -o out.fasta
        seqkit grep -f ${kraken2results} out.fasta -o ${meta.id}_R1.fasta
        rm out.fasta
        seqkit fq2fa ${trimmedreads[1]} -o out.fasta
        seqkit grep -f ${kraken2results} out.fasta -o ${meta.id}_R2.fasta
        rm out.fasta
    fi
    """
}
