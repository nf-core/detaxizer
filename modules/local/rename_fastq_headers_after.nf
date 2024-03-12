process RENAME_FASTQ_HEADERS_AFTER {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::seqkit=2.8.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.8.0--h9ee0642_0':
        'biocontainers/seqkit:2.8.0--h9ee0642_0'}"

    input:
    tuple val(meta), path(fastqfiltered), path(renamedHeaders)

    output:
    tuple val(meta), path('*.fastq.gz'), emit: fastq
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    if [ "$meta.single_end" == "true" ]; then
        gzip -d $renamedHeaders
        seqkit replace -p '^(.+)\$' -r '{kv}' -k *_headers.txt $fastqfiltered -o ${meta.id}_filtered.fastq.gz
        rm *_headers.txt
    else
        gzip -d ${renamedHeaders[0]}
        seqkit replace -p '^(.+)\$' -r '{kv}' -k *_headers_fw.txt ${fastqfiltered[0]} -o ${meta.id}_R1_filtered.fastq.gz
        rm *_headers_fw.txt
        gzip -d ${renamedHeaders[1]}
        seqkit replace -p '^(.+)\$' -r '{kv}' -k *_headers_rv.txt ${fastqfiltered[1]} -o ${meta.id}_R2_filtered.fastq.gz
        rm *_headers_rv.txt
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed -E 's/.*v([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/')
    END_VERSIONS
    """
}
