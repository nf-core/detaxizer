process PREPARE_FASTA4BLASTN {

    container = "https://depot.galaxyproject.org/singularity/seqtk%3A1.4--he4a0461_1"

    input:
        tuple val(meta), path(trimmedreads), path(kraken2results)

    output:
        tuple val(meta), path("*.fasta"), emit: fasta
        tuple val(meta), path("ids*.txt"), emit: ids
        path("versions.yml"), emit: versions

    script:
    """

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: 1.4
    END_VERSIONS

    if [ "$meta.single_end" == "true" ]; then

        awk -F'\t' '{print \$2}' ${kraken2results} > ids.txt
        seqtk subseq ${trimmedreads} ids.txt | seqtk seq -A - > ${meta.id}.fasta    

    else
        awk -F'\t' '{print \$2"/1"}' ${kraken2results} > ids_R1.txt
        awk -F'\t' '{print \$2"/2"}' ${kraken2results} > ids_R2.txt

        seqtk subseq ${trimmedreads[0]} ids_R1.txt | seqtk seq -A - > ${meta.id}_R1.fasta
        seqtk subseq ${trimmedreads[1]} ids_R2.txt | seqtk seq -A - > ${meta.id}_R2.fasta
    fi
    """
}
