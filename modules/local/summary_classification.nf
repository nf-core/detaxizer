process SUMMARY_CLASSIFICATION {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.10.4 pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"
    input:
    tuple val(meta),path(classification)

    output:
    tuple val(meta), path("*.classification_summary.tsv")   , emit: summary
    path("versions.yml")                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python
    import subprocess
    import pandas as pd

    def get_version():
        version_output = subprocess.getoutput('python --version')
        return version_output.split()[1]

    with open("${classification}", 'r') as fp:
        lines = len(fp.readlines())
    if ("${params.classification_kraken2}" == 'false' and "${params.classification_bbduk}" == 'false') or ("${params.classification_kraken2}" == 'true' and "${params.classification_bbduk}" == 'false'):
        classified_dict = {
            "classified with kraken2": [str(lines)]
        }
    elif ("${params.classification_kraken2}" == 'false' and "${params.classification_bbduk}" == 'true'):
        classified_dict = {
            "classified with bbduk": [str(lines)]
        }
    elif ("${params.classification_kraken2}" == 'true' and "${params.classification_bbduk}" == 'true'):
        classified_dict = {
            "classified with kraken2 and bbduk": [str(lines)]
        }
    print(classified_dict)
    df = pd.DataFrame(classified_dict)
    index_name = "${meta.id}".replace("_R1","")
    df.index = [index_name]

    df.to_csv(index_name + ".classification_summary.tsv",sep='\\t')

    # Generate the version.yaml for MultiQC
    with open('versions.yml', 'w') as f:
        f.write(f'"{subprocess.getoutput("echo ${task.process}")}":\\n')
        f.write(f'    python: {get_version()}\\n')
        f.write(f'    pandas: {pd.__version__}\\n')
    """
}
