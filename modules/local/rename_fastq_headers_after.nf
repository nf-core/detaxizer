process RENAME_FASTQ_HEADERS_AFTER {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::python=3.10.4 biopython=1.83 numpy=1.26.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81' :
        'biocontainers/biopython:1.81' }"

    input:
    tuple val(meta), path(fastqfiltered), path(dict)

    output:
    tuple val(meta), path('*.fastq.gz'), emit: fastq
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python
    import Bio
    from Bio import SeqIO, bgzf
    import gzip
    import sys
    import json
    import subprocess

    def get_version():
        version_output = subprocess.getoutput('python --version')
        version = version_output.split()[1]
        return version

    with open("${dict}", 'r') as file:
        headerDict = json.load(file)

    fastq = "${fastqfiltered}".split(" ")
    if len(fastq) == 2:
        with gzip.open(fastq[0], "rt") as handle1, gzip.open(fastq[1], "rt") as handle2:
            with bgzf.BgzfWriter("${meta.id}_R1_filtered.fastq.gz", "wb") as outgz1, bgzf.BgzfWriter("${meta.id}_R2_filtered.fastq.gz", "wb") as outgz2:
                for i,j in zip(SeqIO.parse(handle1, "fastq"),SeqIO.parse(handle2, "fastq")):
                    id_fw = i.id
                    id_rv = j.id
                    if id_fw != id_rv:
                        sys.exit("The IDs did not match. The provided fastq files are either not sorted or a sequence is missing.")
                    lookupHeaders = headerDict[id_fw]
                    i.id = lookupHeaders[0]
                    i.description = ""
                    j.id = lookupHeaders[1]
                    j.description = ""
                    SeqIO.write(sequences=i, handle=outgz1, format="fastq")
                    SeqIO.write(sequences=j, handle=outgz2, format="fastq")

    else:
        with gzip.open(fastq[0], "rt") as handle1:
            with bgzf.BgzfWriter("${meta.id}_filtered.fastq.gz", "wb") as outgz1:
                for i in SeqIO.parse(handle1, "fastq"):
                    id_fw = i.id
                    lookupHeaders = headerDict[id_fw]
                    i.id = lookupHeaders[0]
                    i.description = ""
                    SeqIO.write(sequences=i, handle=outgz1, format="fastq")

    with open('versions.yml', 'w') as f:
        f.write(f'"{subprocess.getoutput("echo ${task.process}")}":\\n')
        f.write(f'    python: {get_version()}\\n')
        f.write(f'    biopython: {Bio.__version__}\\n')
    """
}
