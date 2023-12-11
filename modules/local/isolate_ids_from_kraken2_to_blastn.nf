process ISOLATE_IDS_FROM_KRAKEN2_TO_BLASTN {

        container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10.4' :
        'biocontainers/python:3.10.4' }"
    input:
    tuple val(meta), path(kraken2results), path(tax2filter)


    output:
    tuple val(meta), path('*classified.txt'), emit: classified
    tuple val(meta), path('*ids.txt'), emit: classified_ids
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python
    import subprocess
    tax2filter = []
    with open('${tax2filter}','r') as file:
        for line in file:
            line = line.strip('\\n')
            print(line)
            tax2filter.append(line)

    filterList = []
    with open('${kraken2results}', 'r') as file:
        with open('${meta.id}.classified.txt', 'w') as outfile:
            for line in file:
                line = line.split("\\t")
                for entry in tax2filter:
                    if entry in line[3]:
                        filterList.append(line[1])
                        outfile.write("\\t".join(line))

    with open('${meta.id}.ids.txt', 'w') as outfile:
        for entry in filterList:
            outfile.write(entry+"\\n")

    def get_version():
        version_output = subprocess.getoutput('python --version')
        version = version_output.split()[1]
        return version

    with open('versions.yml', 'w') as f:
        f.write(f'"{subprocess.getoutput("echo ${task.process}")}":\\n')
        f.write(f'    python: {get_version()}\\n')
    """
}
