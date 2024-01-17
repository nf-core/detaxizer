process RENAME_FASTQ_HEADERS_PRE {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::python=3.10.4 biopython=1.83 numpy=1.26.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81' :
        'biocontainers/biopython:1.81' }"
    input:
    tuple val(meta), path(inputfastq)

    output:
    tuple val(meta), path('*_headers.json'), emit: headers
    tuple val(meta), path('*.fastq.gz'), emit: fastq
    path "versions.yml", emit: versions

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


    def renameReadsPaired(reads):
        read_fw = reads[0]
        read_rv = reads[1]
        read_dict = {}
        if "/1" in read_fw and "/2" in read_rv and " " not in read_fw and " " not in read_rv:
            read_fw_stripped = read_fw.strip("1").strip("/")
            read_rv_stripped = read_rv.strip("2").strip("/")
            if read_fw_stripped != read_rv_stripped:
                sys.exit("Read IDs were not matching! Please provide matching headers.")
            else:
                read_dict[read_fw_stripped] = [read_fw, read_rv]
                read_fw_stripped = read_fw_stripped + " 1:N:10:"
                read_rv_stripped = read_rv_stripped + " 2:N:10:"
                read_renamed = [read_fw_stripped,read_rv_stripped]
        elif "/1" in read_fw and "/2" in read_rv and " " in read_fw and " " in read_rv:
            read_fw_stripped = read_fw.strip("1").strip("/").split(" ")[0]
            read_rv_stripped = read_rv.strip("2").strip("/").split(" ")[0]
            if read_fw_stripped != read_rv_stripped:
                sys.exit("Read IDs were not matching! Please provide matching headers.")
            else:
                read_dict[read_fw_stripped] = [read_fw, read_rv]
                read_fw_stripped = read_fw_stripped + " 1:N:10:"
                read_rv_stripped = read_rv_stripped + " 2:N:10:"
                read_renamed = [read_fw_stripped,read_rv_stripped]
        elif " 1:" in read_fw and " 2:" in read_rv:
            read_fw_stripped = read_fw.split(" ")[0]
            read_rv_stripped = read_rv.split(" ")[0]
            if read_fw_stripped != read_rv_stripped:
                sys.exit("Read IDs were not matching! Please provide matching headers.")
            else:
                read_dict[read_fw_stripped] = [read_fw, read_rv]
                read_fw_stripped = read_fw_stripped + " 1:N:10:"
                read_rv_stripped = read_rv_stripped + " 2:N:10:"
                read_renamed = [read_fw_stripped,read_rv_stripped]
        else:
            sys.exit("The headers were not matching the patterns!")
        return (read_dict,read_renamed)

    def renameReadSingle(read):
        read_dict = {}
        if "/1" in read and " " not in read:
            read_fw_stripped = read.strip("1").strip("/")
            read_dict[read_fw_stripped] = [read]
            read_fw_stripped = read_fw_stripped + " 1:N:10:"
            read_renamed = [read_fw_stripped]
        elif "/1" in read and " " in read:
            read_fw_stripped = read.strip("1").strip("/").split(" ")[0]
            read_dict[read_fw_stripped] = [read]
            read_fw_stripped = read_fw_stripped + " 1:N:10:"
            read_renamed = [read_fw_stripped]
        elif " 1:" in read:
            read_fw_stripped = read.split(" ")[0]
            read_dict[read_fw_stripped] = [read]
            read_fw_stripped = read_fw_stripped + " 1:N:10:"
            read_renamed = [read_fw_stripped]
        elif " " in read:
            read_fw_stripped = read.split(" ")[0]
            read_dict[read_fw_stripped] = [read]
            read_fw_stripped = read_fw_stripped + " 1:N:10:"
            read_renamed = [read_fw_stripped]
        else:
            sys.exit("The headers were not matching the patterns!")
        return (read_dict,read_renamed)


    fastq = "${inputfastq}".split(" ")
    if len(fastq) == 2:
        renamed = {}
        with gzip.open(fastq[0], "rt") as handle1, gzip.open(fastq[1], "rt") as handle2:
            with bgzf.BgzfWriter("${meta.id}_R1_renamed.fastq.gz", "wb") as outgz1, bgzf.BgzfWriter("${meta.id}_R2_renamed.fastq.gz", "wb") as outgz2:
                for i,j in zip(SeqIO.parse(handle1, "fastq"),SeqIO.parse(handle2, "fastq")):
                    header_fw = i.description
                    header_rv = j.description
                    headers = (header_fw,header_rv)
                    headers = renameReadsPaired(headers)
                    renamed.update(headers[0])
                    i.description = headers[1][0]
                    i.id = headers[1][0].split(" ")[0]
                    j.description = headers[1][1]
                    j.id = headers[1][1].split(" ")[0]
                    SeqIO.write(sequences=i, handle=outgz1, format="fastq")
                    SeqIO.write(sequences=j, handle=outgz2, format="fastq")
        with open("${meta.id}_headers.json", "w") as outfile:
            json.dump(renamed, outfile)
    else:
        renamed = {}
        with gzip.open(fastq[0], "rt") as handle1:
            with bgzf.BgzfWriter("${meta.id}_renamed.fastq.gz", "wb") as outgz1:
                for i in SeqIO.parse(handle1, "fastq"):
                    header_fw = i.description
                    headers = renameReadSingle(header_fw)
                    renamed.update(headers[0])
                    i.description = headers[1][0]
                    i.id = headers[1][0].split(" ")[0]
                    SeqIO.write(sequences=i, handle=outgz1, format="fastq")
        with open("${meta.id}_headers.json", "w") as outfile:
            json.dump(renamed, outfile)

    def get_version():
        version_output = subprocess.getoutput('python --version')
        version = version_output.split()[1]
        return version

    with open('versions.yml', 'w') as f:
        f.write(f'"{subprocess.getoutput("echo ${task.process}")}":\\n')
        f.write(f'    python: {get_version()}\\n')
        f.write(f'    biopython: {Bio.__version__}\\n')
    """
}

process RENAME_FASTQ_HEADERS_AFTER {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::biopython=1.83"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81' :
        'biocontainers/biopython:1.81' }"

    input:
    tuple val(meta), path(fastqfiltered), path(dict)

    output:
    tuple val(meta), path('*.fastq.gz'), emit: fastq
    path "versions.yml", emit: versions

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

    def get_version():
        version_output = subprocess.getoutput('python --version')
        version = version_output.split()[1]
        return version

    with open('versions.yml', 'w') as f:
        f.write(f'"{subprocess.getoutput("echo ${task.process}")}":\\n')
        f.write(f'    python: {get_version()}\\n')
        f.write(f'    biopython: {Bio.__version__}\\n')
    """
}
