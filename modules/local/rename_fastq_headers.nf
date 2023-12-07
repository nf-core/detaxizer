process RENAME_FASTQ_HEADERS {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81' :
        'biocontainers/biopython:1.81' }"
    input:
    tuple val(meta), path(inputfastq)

    output:
    tuple val(meta), path('*_headers.json'), emit: headers
    tuple val(meta), path('*.fastq.gz'), emit: fastq
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python
    from Bio import SeqIO, bgzf
    import gzip
    import sys
    import json
    import re

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
            read_dict[read_fw_stripped] = [read, read_rv]
            read_dict[read_fw_stripped] = [read]
            read_fw_stripped = read_fw_stripped + " 1:N:10:"
            read_renamed = [read_fw_stripped]
        elif " 1:" in read_fw and " 2:" in read_rv:
            read_fw_stripped = read.split(" ")[0]
            read_dict[read_fw_stripped] = [read]
            read_fw_stripped = read_fw_stripped + " 1:N:10:"
            read_renamed = [read_fw_stripped]
        else:
            sys.exit("The headers were not matching the patterns!")
        return (read_dict,read_renamed)


    fastq = "${inputfastq}".split(" ")
    if len(fastq) == 2:
        with gzip.open(fastq[0], "rt") as handle:
            reads_fw = list(SeqIO.parse(handle, "fastq"))
        with gzip.open(fastq[1], "rt") as handle:
            reads_rv = list(SeqIO.parse(handle, "fastq"))
        reads = zip(reads_fw,reads_rv)
        renamed = {}
        new_reads = []
        for read_fw_rv in reads:
            header_fw = read_fw_rv[0].description
            header_rv = read_fw_rv[1].description
            headers = (header_fw,header_rv)
            headers = renameReadsPaired(headers)
            renamed.update(headers[0])
            read_fw_rv[0].description = headers[1][0]
            read_fw_rv[0].id = headers[1][0].split(" ")[0]
            read_fw_rv[1].description = headers[1][1]
            read_fw_rv[1].id = headers[1][1].split(" ")[0]
            new_reads.append(read_fw_rv)
        del reads
        with open("${meta.id}_headers.json", "w") as outfile:
            json.dump(renamed, outfile)
        del renamed
        reads_fw_renamed,reads_rv_renamed = zip(*new_reads)
        with bgzf.BgzfWriter("${meta.id}_R1_renamed.fastq.gz", "wb") as outgz:
            SeqIO.write(sequences=reads_fw_renamed, handle=outgz, format="fastq")
        del reads_fw_renamed
        with bgzf.BgzfWriter("${meta.id}_R2_renamed.fastq.gz", "wb") as outgz:
            SeqIO.write(sequences=reads_rv_renamed, handle=outgz, format="fastq")
        del reads_rv_renamed
    else:
        with gzip.open(fastq[0], "rt") as handle:
            reads_fw = list(SeqIO.parse(handle, "fastq"))
        renamed = {}
        new_reads = []
        for read_fw in reads_fw:
            header_fw = read_fw.description
            headers = renameReadSingle(header_fw)
            renamed.update(headers[0])
            read_fw.description = headers[1][0]
            read_fw.id = headers[1][0].split(" ")[0]
            new_reads.append(read_fw)
        del reads_fw
        with open("${meta.id}_headers.json", "w") as outfile:
            json.dump(renamed, outfile)
        del renamed
        with bgzf.BgzfWriter("${meta.id}_renamed.fastq.gz", "wb") as outgz:
            SeqIO.write(sequences=new_reads, handle=outgz, format="fastq")
        del new_reads
    import subprocess
    def get_version():
        version_output = subprocess.getoutput('python --version')
        version = version_output.split()[1]
        return version

    with open('versions.yml', 'w') as f:
        f.write(f'"{subprocess.getoutput("echo ${task.process}")}":\\n')
        f.write(f'    python: {get_version()}\\n')
    """
}
